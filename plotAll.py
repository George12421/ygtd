# 综合对比版: C-Path vs PageRank vs NSG (含时间统计表格图)
# 适配多场景目录结构 (scenario_L*_R*)
# 支持自由开关算法、多进程加速、统一测试集并在同一张图上绘制对比曲线

import pandas as pd
import numpy as np
import math
import os
import scipy.io as sio
import glob
import time
import heapq
import networkx as nx
from collections import deque, defaultdict
import matplotlib.pyplot as plt
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import sys

# =========================================================================
# 0. 配置区 (算法开关与参数)
# =========================================================================

# 自由开关需要参与对比的算法
ENABLE_ALGORITHMS = {
    'C-Path': True,  # 原论文 Baseline
    'PageRank': True,  # 改进方案一
    'NSG': True  # 改进方案二
}

# 测试的 Top-K 比例列表
K_PERCENT_LIST = [0.10, 0.20, 0.30, 0.40, 0.50]
VISUALIZATION_THRESHOLD = 0.05
MAX_VALIDATION_STAGE = 3             #14-7  39-5  118-3
CASE_NAME = 'IEEE118'

# NSG 专用参数
NSG_MAX_DEGREE = 1000  # 稀疏化限制
NSG_EXPANSION = 1000  # 最大搜索扩展深度

# 新增：掩码值控制变量列表
MASK_VALUES = [0.1, 0.3, 0.5, 0.7, 1.0]

# 路径配置 (注意 BASE_PATH_GT 指向多场景测试集的根目录)
BASE_PATH_B_MATRIX = Path(fr'./data/dc/{CASE_NAME}/causal_results_csv')
BASE_PATH_GT = fr'./data/dc/{CASE_NAME}/test_cascade_data'  # 修改为你的多场景文件夹名
PLOT_OUTPUT_DIR = Path(fr'./data/dc/{CASE_NAME}')
PLOT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# =========================================================================
# 1. 算法核心库 (集成了三种算法)
# =========================================================================

# --- 算法A: C-Path (BFS) ---
def find_all_paths_bfs(adj_matrix, start, end, max_len=3):
    paths = []
    queue = deque([(start, [start])])
    while queue:
        curr_node, path = queue.popleft()
        if len(path) - 1 >= max_len: continue
        neighbors = np.where(adj_matrix[:, curr_node] != 0)[0]
        for neighbor in neighbors:
            if neighbor == end:
                paths.append(path + [neighbor])
            elif neighbor not in path and len(path) < max_len + 1:
                queue.append((neighbor, path + [neighbor]))
    return paths


def calculate_total_effect_dm(B_modified, start_node, end_node, max_len=3):
    all_paths = find_all_paths_bfs(B_modified, start_node, end_node, max_len)
    if not all_paths: return 0.0
    total_effect = 0.0
    for path in all_paths:
        path_score = 1.0
        for i in range(len(path) - 1): path_score *= B_modified[path[i + 1], path[i]]
        total_effect += path_score
    return total_effect


def predict_next_cascade_cpath(B_modified, most_recent_fault, all_failed_lines, k_percent):
    num_lines = B_modified.shape[0]
    candidates = set(range(num_lines)) - all_failed_lines
    if not candidates: return set()
    candidate_scores = {j: calculate_total_effect_dm(B_modified, most_recent_fault, j) for j in candidates}
    total_abs_sum = sum(abs(score) for score in candidate_scores.values())
    final_scores = {j: (abs(score) / total_abs_sum) if total_abs_sum > 0 else 0.0 for j, score in
                    candidate_scores.items()}
    num_to_select = int(math.ceil(k_percent * num_lines))
    return {j for j, _ in sorted(final_scores.items(), key=lambda item: item[1], reverse=True)[:num_to_select]}


# --- 算法B: PageRank ---
def predict_next_cascade_pagerank(B_modified, most_recent_fault, all_failed_lines, k_percent, damping=0.85):
    num_lines = B_modified.shape[0]
    G = nx.from_numpy_array(np.abs(B_modified).T, create_using=nx.DiGraph)
    personalization = dict.fromkeys(G.nodes, 0.0)
    if most_recent_fault in G:
        personalization[most_recent_fault] = 1.0
    else:
        return set()
    try:
        scores = nx.pagerank(G, alpha=damping, personalization=personalization, weight='weight')
    except:
        return set()
    candidates = set(range(num_lines)) - all_failed_lines
    final_scores = {node: score for node, score in scores.items() if node in candidates}
    num_to_select = int(math.ceil(k_percent * num_lines))
    return {j for j, _ in sorted(final_scores.items(), key=lambda item: item[1], reverse=True)[:num_to_select]}


# --- 算法C: NSG ---
def predict_next_cascade_nsg(B_modified, most_recent_fault, all_failed_lines, k_percent, max_degree=20,
                             expansion_depth=200):
    num_lines = B_modified.shape[0]
    A = np.abs(B_modified).T

    def get_top_k_neighbors(node_idx, k):
        if node_idx >= num_lines: return []
        row = A[node_idx]
        non_zero_indices = np.nonzero(row)[0]
        if len(non_zero_indices) == 0: return []
        weights = row[non_zero_indices]
        if len(non_zero_indices) <= k: return [(idx, w) for idx, w in zip(non_zero_indices, weights)]
        top_k_idx_local = np.argpartition(weights, -k)[-k:]
        return [(non_zero_indices[i], weights[i]) for i in top_k_idx_local]

    entry_point = most_recent_fault
    visited, candidate_pool, final_scores = {entry_point}, [], defaultdict(float)
    for neighbor, weight in get_top_k_neighbors(entry_point, max_degree):
        if neighbor not in visited and neighbor not in all_failed_lines:
            heapq.heappush(candidate_pool, (-weight, neighbor))
            visited.add(neighbor)

    expanded_nodes = 0
    while candidate_pool and expanded_nodes < expansion_depth:
        neg_weight, current_node = heapq.heappop(candidate_pool)
        current_weight = -neg_weight
        final_scores[current_node] += current_weight
        expanded_nodes += 1
        for neighbor, edge_weight in get_top_k_neighbors(current_node, max_degree):
            if neighbor not in visited and neighbor not in all_failed_lines:
                heapq.heappush(candidate_pool, (-(current_weight * edge_weight), neighbor))
                visited.add(neighbor)

    candidates = set(range(num_lines)) - all_failed_lines
    for j in candidates:
        if j not in final_scores: final_scores[j] = 0.0
    num_to_select = int(math.ceil(k_percent * num_lines))
    return {j for j, _ in sorted(final_scores.items(), key=lambda item: item[1], reverse=True)[:num_to_select]}


# =========================================================================
# 2. 通用工具与验证逻辑 (多场景适配)
# =========================================================================

def load_b_matrix(line_index, base_path, num_lines):
    found_files = glob.glob(os.path.join(base_path, f"causal_B_anomaly_F{line_index:04d}_R*.csv"))
    filepath = found_files[0] if found_files else os.path.join(base_path, f"causal_B_anomaly_F{line_index:04d}.csv")
    return pd.read_csv(filepath).to_numpy()


def get_modified_b_matrix(fault_chain, base_path, num_lines, threshold):
    most_recent_fault = fault_chain[-1]
    B_modified = load_b_matrix(most_recent_fault, base_path, num_lines).copy()
    for i in set(fault_chain[:-1]):
        if 0 <= i < num_lines: B_modified[i, :] = 0
    B_modified[np.abs(B_modified) < threshold] = 0.0
    return B_modified, most_recent_fault, set(fault_chain)


def parse_ground_truth_chain(mat_filepath, initial_fault_index):
    active_lines_data = sio.loadmat(mat_filepath)['activeLinesData']
    num_stages_M = active_lines_data.shape[1]
    if num_stages_M < 2: return [initial_fault_index]
    active_set = [set(np.where(active_lines_data[:, m] == 1)[0]) for m in range(num_stages_M)]
    gt_chain = [initial_fault_index]
    for m in range(2, num_stages_M):
        failed_lines = active_set[m - 1] - active_set[m]
        if failed_lines:
            gt_chain.append(list(failed_lines)[0])
        else:
            break
    return gt_chain


def validate_single_sequence(gt_chain, base_b_matrix_path, num_lines, k_percent, threshold, max_stage, algo_name):
    if len(gt_chain) < 2: return np.nan, []
    total_pred = correct_pred = 0
    exec_times = []

    for m in range(1, len(gt_chain)):
        if m + 1 > max_stage: break
        try:
            B_mod, most_recent, all_failed = get_modified_b_matrix(gt_chain[:m], base_b_matrix_path, num_lines,
                                                                   threshold)

            start_time = time.time()
            if algo_name == 'C-Path':
                predicted_set = predict_next_cascade_cpath(B_mod, most_recent, all_failed, k_percent)
            elif algo_name == 'PageRank':
                predicted_set = predict_next_cascade_pagerank(B_mod, most_recent, all_failed, k_percent)
            elif algo_name == 'NSG':
                predicted_set = predict_next_cascade_nsg(B_mod, most_recent, all_failed, k_percent,
                                                         max_degree=NSG_MAX_DEGREE, expansion_depth=NSG_EXPANSION)

            exec_times.append(time.time() - start_time)

            if gt_chain[m] in predicted_set: correct_pred += 1
            total_pred += 1
        except:
            break

    precision = correct_pred / total_pred if total_pred > 0 else np.nan
    return precision, exec_times


def process_single_file_wrapper(args):
    """ 多进程包装函数 (已适配多场景的 3 元组解包) """
    scenario_name, idx, gt_file, k_val, base_b, n_lines, thresh, mx_stage, algo = args
    return validate_single_sequence(parse_ground_truth_chain(gt_file, idx), base_b, n_lines, k_val, thresh, mx_stage,
                                    algo)


def discover_gt_files(base_gt_path):
    """ 穿透扫描多场景目录 (scenario_L*_R*) 下的所有序列 """
    files = []
    if not os.path.exists(base_gt_path):
        return files

    for scenario_dir in os.listdir(base_gt_path):
        scenario_path = os.path.join(base_gt_path, scenario_dir)
        if not os.path.isdir(scenario_path) or not scenario_dir.startswith('scenario_'):
            continue

        for dir_name in os.listdir(scenario_path):
            full_dir_path = os.path.join(scenario_path, dir_name)
            if os.path.isdir(full_dir_path) and dir_name.startswith('failure_'):
                try:
                    idx = int(dir_name[8:])
                    for mat_path in glob.glob(os.path.join(full_dir_path, 'cascade_seq_*.mat')):
                        files.append((scenario_dir, idx, mat_path))
                except ValueError:
                    pass
    # 按场景名和初始故障排序
    return sorted(files, key=lambda x: (x[0], x[1]))


# =========================================================================
# 3. 主程序执行与多曲线、表格绘图
# =========================================================================
if __name__ == "__main__":

    # 规模探测
    NUM_LINES = None
    for mask in MASK_VALUES:
        mask_dir = f'mask_{mask}'
        probe_path = BASE_PATH_B_MATRIX / mask_dir
        probe = glob.glob(os.path.join(probe_path, 'causal_B_anomaly_F*.csv'))
        if probe:
            NUM_LINES = pd.read_csv(probe[0]).shape[0]
            break
    if NUM_LINES is None:
        print(f"错误: 在所有掩码文件夹下未找到 B 矩阵")
        exit(1)

    files_to_validate = discover_gt_files(BASE_PATH_GT)
    if not files_to_validate:
        print(f"错误: 在 {BASE_PATH_GT} 下未找到任何多场景级联序列文件")
        exit(1)

    print(f"成功扫描到 {len(files_to_validate)} 个跨场景测试序列。")

    max_workers = max(1, multiprocessing.cpu_count() - 4)

    for mask in MASK_VALUES:
        mask_dir = f'mask_{mask}'
        base_b_path = BASE_PATH_B_MATRIX / mask_dir
        print(f"\n{'=' * 60}")
        print(f" 开始处理掩码值: {mask} (K={min(K_PERCENT_LIST) * 100:.0f}% ~ {max(K_PERCENT_LIST) * 100:.0f}%)")
        print(f"{'=' * 60}")

        all_precisions = {}
        all_time_metrics = {}

        for algo_name, is_enabled in ENABLE_ALGORITHMS.items():
            if not is_enabled: continue

            print(f"\n[启动] 正在测试算法: {algo_name}")
            algo_precisions = []
            algo_all_single_step_times = []
            algo_start_time = time.time()

            for k_val in K_PERCENT_LIST:
                # 组装多进程任务参数
                task_args = [
                    (scenario, idx, gt, k_val, str(base_b_path), NUM_LINES, VISUALIZATION_THRESHOLD, MAX_VALIDATION_STAGE,
                     algo_name)
                    for scenario, idx, gt in files_to_validate]

                valid_scores = []
                completed = 0

                with ProcessPoolExecutor(max_workers=max_workers) as executor:
                    for future in as_completed([executor.submit(process_single_file_wrapper, arg) for arg in task_args]):
                        completed += 1
                        sys.stdout.write(
                            f"\r  --> K={k_val * 100:2.0f}% | 进度: {completed}/{len(task_args)} ({completed / len(task_args) * 100:.1f}%)")
                        sys.stdout.flush()

                        score, single_step_times = future.result()
                        if not np.isnan(score):
                            valid_scores.append(score)
                            algo_all_single_step_times.extend(single_step_times)

                avg_p = np.mean(valid_scores) if valid_scores else 0.0
                algo_precisions.append(avg_p)
                print(f" | 平均精确度: {avg_p * 100:.2f}%")

            algo_total_time = time.time() - algo_start_time
            avg_single_step_time = np.mean(algo_all_single_step_times) if algo_all_single_step_times else 0.0

            all_precisions[algo_name] = algo_precisions
            all_time_metrics[algo_name] = {
                'avg_single_step': avg_single_step_time,
                'total_exec_time': algo_total_time
            }

        # =========================================================================
        # 4. 绘图逻辑 1: 精确度对比柱状图
        # =========================================================================
        print("\n============================================================")
        print(f" 开始绘制掩码值 {mask} 的联合对比柱状图与时间统计表格...")

        plt.rcParams.update({'font.family': 'serif', 'font.size': 14})

        plt.figure(figsize=(10, 7))
        x = np.arange(len(K_PERCENT_LIST))
        width = 0.25  # 柱宽

        styles = {
            'C-Path': {'color': '#2ca02c', 'label': 'C-Path (Baseline)'},
            'PageRank': {'color': '#C23531', 'label': 'Personalized PageRank'},
            'NSG': {'color': '#2F4554', 'label': 'NSG Greedy Search'}
        }

        for i, algo_name in enumerate(all_precisions):
            y_values = [p * 100 for p in all_precisions[algo_name]]
            plt.bar(x + i * width, y_values, width, label=styles[algo_name]['label'], color=styles[algo_name]['color'])

            # 添加数值标签
            for j, y in enumerate(y_values):
                plt.text(x[j] + i * width, y + 0.5, f'{y:.2f}', ha='center', va='bottom', fontsize=11, color=styles[algo_name]['color'])

        plt.xlabel('Top-K Selection ($\kappa$ %)', fontsize=16, fontweight='bold')
        plt.ylabel('Overall Prediction Precision (%)', fontsize=16, fontweight='bold')
        plt.title(f'Algorithm Comparison on Precision ({CASE_NAME}, Mask={mask})', fontsize=18, fontweight='bold')

        plt.xticks(x + width, [int(k * 100) for k in K_PERCENT_LIST])
        all_y = [p * 100 for precisions in all_precisions.values() for p in precisions]
        plt.ylim(max(0, min(all_y) - 5), min(100, max(all_y) + 10))
        plt.grid(True, linestyle='--', alpha=0.6, axis='y')
        plt.legend(loc='upper left', fontsize=12, frameon=True, shadow=True)

        algo_str = "_".join(all_precisions.keys())
        plot_filepath_bar = PLOT_OUTPUT_DIR / f'precision_comparison_{CASE_NAME}_mask_{mask}_multiscenario_{algo_str}.png'
        plt.tight_layout()
        plt.savefig(plot_filepath_bar, dpi=300, bbox_inches='tight')
        plt.show()

        # =========================================================================
        # 5. 绘图逻辑 2: 时间统计表格图
        # =========================================================================
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.axis('tight')
        ax.axis('off')

        table_data = []
        row_labels = []

        for algo_name, metrics in all_time_metrics.items():
            row_labels.append(styles[algo_name]['label'])
            table_data.append([f"{metrics['avg_single_step']:.5f} s", f"{metrics['total_exec_time']:.2f} s"])

        col_labels = ['Avg. Single Step Time', 'Total Execution Time']

        table = ax.table(cellText=table_data,
                         rowLabels=row_labels,
                         colLabels=col_labels,
                         loc='center',
                         cellLoc='center')

        table.auto_set_font_size(False)
        table.set_fontsize(14)
        table.scale(1, 2.5)

        for (row, col), cell in table.get_celld().items():
            if row == 0 or col == -1:
                cell.set_text_props(weight='bold')

        plt.title(f"Computational Efficiency Comparison ({CASE_NAME}, Mask={mask})", fontsize=16, fontweight='bold', pad=20)

        plot_filepath_table = PLOT_OUTPUT_DIR / f'time_table_comparison_{CASE_NAME}_mask_{mask}_multiscenario_{algo_str}.png'
        plt.tight_layout()
        plt.savefig(plot_filepath_table, dpi=300, bbox_inches='tight')
        plt.show()

        print(f"\n掩码值 {mask} 绘图全部完成！")
        print(f" 1. 联合精确度对比柱状图已保存至: {plot_filepath_bar}")
        print(f" 2. 时间统计表格图已保存至      : {plot_filepath_table}")
