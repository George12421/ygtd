import pandapower as pp
import pandapower.networks as pn
import numpy as np
import pandas as pd
from scipy.linalg import eigh, inv
from scipy.optimize import linear_sum_assignment
import sys
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
from pathlib import Path
import warnings

# 字体设置
try:
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman'] + plt.rcParams['font.serif'],
        'mathtext.fontset': 'stix'
    })
    print("字体 'Times New Roman' 设置成功。")
except:
    print("警告： 'Times New Roman' 字体未找到，将使用默认衬线字体。")


# --- 1. 配置区 ---

CASE_NAME = 'IEEE118'
ROOT_PATH = Path('D:/xuexi\dianwangjilian/yinguoyuce/ieee_2')  # 请根据实际情况修改
BASE_DIR = ROOT_PATH / f'data/dc/{CASE_NAME}'
NPZ_INPUT_DIR = BASE_DIR / 'anomaly_data_npz'
CSV_OUTPUT_DIR = BASE_DIR / 'causal_results_csv'
PLOT_OUTPUT_DIR = BASE_DIR / 'causal_results_png'

# 新增：掩码值控制变量列表
MASK_VALUES = [0.1, 0.3, 0.5, 0.7, 1.0]

# (!!) NEW: 指定物理拓扑文件的路径
# 注意：这里的拓扑必须是“线路-线路”的邻接矩阵 (Line Connectivity Matrix)
# 也就是：如果两条线路共享同一个母线，则它们是相连的 (值为1)，否则为0。
# 维度必须是 (N_lines, N_lines)
TOPOLOGY_FILE = BASE_DIR / 'physical_line_topology.csv'


def generate_edge_adjacency(output_path, case_name='IEEE118'):
    """
    通用版：根据 case_name 生成线路物理邻接矩阵。
    【修改版】：自动合并连接相同母线对的平行线路 (Parallel Lines)。
    """
    print(f"正在从 pandapower 生成 {case_name} 拓扑矩阵 (含平行线路合并)...")

    # 1. 根据配置加载对应的案例
    try:
        if case_name == 'IEEE14':
            net = pn.case14()
        elif case_name == 'IEEE39':
            net = pn.case39()
        elif case_name == 'IEEE89':
            net = pn.case89pegase()
        elif case_name == 'IEEE118':
            net = pn.case118()
        else:
            print(f"Error: 未知的 CASE_NAME: {case_name}")
            return None
    except Exception as e:
        print(f"Error: 无法加载 pandapower 案例。错误: {e}")
        return None

    # 2. 提取所有分支 (Branches)
    lines = net.line[['from_bus', 'to_bus']].copy()
    trafos = net.trafo[['hv_bus', 'lv_bus']].copy()
    trafos.columns = ['from_bus', 'to_bus']

    # 原始所有边
    raw_branches = pd.concat([lines, trafos], ignore_index=True)
    print(f"  - 原始分支总数 (含平行线): {len(raw_branches)}")

    # 3. 【核心修改】合并平行线路
    # 逻辑：如果 Line A 连接 (1, 2)，Line B 也连接 (1, 2)，则视为同一条边。

    # 3.1 标准化方向：确保 (from, to) 总是 (小, 大)，忽略电流方向差异
    # 这样 (Bus1, Bus2) 和 (Bus2, Bus1) 都会变成 (Bus1, Bus2)
    p_from = raw_branches[['from_bus', 'to_bus']].min(axis=1)
    p_to = raw_branches[['from_bus', 'to_bus']].max(axis=1)

    # 3.2 创建标准化后的临时 DataFrame
    normalized_df = pd.DataFrame({'u': p_from, 'v': p_to})

    # 3.3 去重：只保留唯一的母线对
    # keep='first' 保持第一次出现的顺序，这对索引对齐很重要
    merged_branches = normalized_df.drop_duplicates(keep='first').reset_index(drop=True)

    n_merged = len(merged_branches)
    print(f"  - 合并平行线路后的有效分支数: {n_merged}")
    print(f"    (已移除 {len(raw_branches) - n_merged} 条平行冗余边)")

    # 4. 构建邻接矩阵 (基于合并后的边)
    adj_matrix = np.zeros((n_merged, n_merged), dtype=int)

    # 将每一行转为集合，方便判断交集 (即是否共享母线)
    branch_sets = [set(row) for row in merged_branches.itertuples(index=False)]

    for i in range(n_merged):
        for j in range(i + 1, n_merged):
            # 如果两条边的端点集合有交集（即共享至少一个 Bus）
            if not branch_sets[i].isdisjoint(branch_sets[j]):
                adj_matrix[i, j] = 1
                adj_matrix[j, i] = 1

    # 5. 保存
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(adj_matrix).to_csv(output_path, header=False, index=False)
    print(f"  - 物理拓扑矩阵已保存至: {output_path}")

    return adj_matrix

# --- 2. SparseICA 辅助函数 (保持不变) ---

def whiten_data(X):
    X_mean = X.mean(axis=1, keepdims=True)
    X_centered = X - X_mean
    cov = np.cov(X_centered)
    d, V = eigh(cov)
    D_inv_sqrt = np.diag(1.0 / np.sqrt(d + 1e-6))
    whitening_matrix = D_inv_sqrt @ V.T
    X_whitened = whitening_matrix @ X_centered
    return X_whitened, whitening_matrix


def gscad_derivative(W, lambda_, lambda1, a, m):
    beta = np.abs(W)
    term1 = np.tanh(m * beta) * (beta <= lambda1)
    term2_numerator = np.maximum(0, a * lambda1 - beta)
    term2_denominator = (a - 1) * lambda1
    if term2_denominator == 0: term2_denominator = 1e-6
    term2 = np.tanh(m * lambda1) * (term2_numerator / term2_denominator) * (beta > lambda1)
    p_prime_abs = lambda_ * (term1 + term2)
    return p_prime_abs * np.sign(W)


def score_function_tanh(Y):
    return np.tanh(Y)


def sparse_ica(X, lambda_=0.05, lambda1=0.075, a=3.7, m=200,
               n_iter=500, learning_rate=0.01, threshold=0.005,
               tolerance=1e-6):
    n_features, n_samples = X.shape
    I = np.eye(n_features)
    X_whitened, whitening_matrix = whiten_data(X)
    W = np.eye(n_features)
    converged = False
    for i in range(n_iter):
        W_old = W.copy()
        Y = W @ X_whitened
        Psi_Y = score_function_tanh(Y)
        E_psi_y_yT = (Psi_Y @ Y.T) / n_samples
        P_prime_W = gscad_derivative(W, lambda_, lambda1, a, m)
        gradient = I - E_psi_y_yT - (P_prime_W @ W.T)
        delta_W = gradient @ W
        W = W + learning_rate * delta_W
        try:
            U, S, Vh = np.linalg.svd(W @ W.T)
            W = (U @ np.diag(1.0 / np.sqrt(S)) @ U.T) @ W
        except np.linalg.LinAlgError:
            pass
        change = np.max(np.abs(W - W_old))
        if change < tolerance:
            # print(f"    - SparseICA 在第 {i + 1} 次迭代收敛")
            converged = True
            break
    W_ica = W @ whitening_matrix
    W_sparse = W_ica.copy()
    W_sparse[np.abs(W_sparse) < threshold] = 0.0
    return W_sparse


# --- 3. (!!) NEW: 物理先验与拓扑加载函数 ---

def load_line_topology(csv_path, expected_dim):
    """
    加载物理拓扑矩阵 (线路-线路连接矩阵)。
    如果文件不存在，返回 None。
    """
    path = Path(csv_path)
    if not path.exists():
        print(f"警告: 物理拓扑文件未找到: {path}")
        print("      将退化为纯数据驱动模式 (不使用物理先验)。")
        return None

    try:
        # 假设没有表头，如果你的csv有表头，请去掉 header=None
        adj = pd.read_csv(path, header=None).values
        if adj.shape != (expected_dim, expected_dim):
            print(f"警告: 物理拓扑维度 {adj.shape} 与数据维度 ({expected_dim}, {expected_dim}) 不匹配。")
            return None

        # 确保是 0/1 矩阵
        adj[adj != 0] = 1
        return adj
    except Exception as e:
        print(f"错误: 读取物理拓扑失败: {e}")
        return None


def apply_physics_prior(B_data, adj_phys, penalty_factor=0.1):  # [cite: 7]
    """
    实施方案一：物理引导的软掩膜 (Physics-Guided Soft Masking)。

    参数:
        B_data: 数据驱动算出的因果矩阵 (N x N)
        adj_phys: 物理邻接矩阵 (N x N), 1表示相连, 0表示不相连
        penalty_factor: 非物理连接边的保留比例 (0.0-1.0)。
                        0.0 表示强制切断 (Hard Constraint)，
                        0.3 表示衰减 70% (Soft Constraint)，
                        1.0 表示不使用先验。
    返回:
        B_refined: 修正后的因果矩阵
    """
    if adj_phys is None:
        return B_data

    # 创建掩膜:
    # 如果 adj_phys 为 1 (物理相连)，mask 为 1.0 (完全保留)
    # 如果 adj_phys 为 0 (物理断开)，mask 为 penalty_factor (抑制)
    mask = np.ones_like(B_data)
    mask[adj_phys == 0] = penalty_factor

    # 应用掩膜
    B_refined = B_data * mask

    return B_refined


# --- 4. 修改后的分析函数 ---

def analyze_causal_matrix(npz_file_path, csv_dir, plot_dir, physics_adj=None, penalty_factor=0.3):
    """
    为单个 .npz 文件运行 ICA 和因果分析 (集成物理先验)。
    完美复现 Cyclic-LiNGAM 论文中的 B 矩阵求解逻辑。
    """
    base_name = npz_file_path.stem
    csv_output_path = csv_dir / f"causal_B_{base_name}.csv"
    plot_output_path = plot_dir / f"causal_graph_{base_name}.png"

    print(f"\n--- Processing: {npz_file_path.name} ---")

    try:
        data = np.load(str(npz_file_path))
        if 'anomalyMatrix' not in data.files:
            return

        s_train_original = data['anomalyMatrix'].copy()
        N_lines, N_loads = s_train_original.shape
        N = N_lines

        # --- Step 1: SparseICA ---
        print("  Running SparseICA (Step 1)...")
        W_ica = sparse_ica(s_train_original, threshold=0.02)

        # --- Step 2: 最佳分配 (Finding the best Assignment) ---
        # 严格遵守论文公式 (14): 最小化对角线绝对值的倒数之和
        print("  Finding best assignment (Step 2)...")
        epsilon_val = 1e-9
        # 构建代价矩阵：取绝对值的倒数
        cost_matrix = 1.0 / (np.abs(W_ica) + epsilon_val)

        # 使用匈牙利算法求解最小代价指派
        row_ind, col_ind = linear_sum_assignment(cost_matrix)

        # 行置换逻辑：
        # col_ind[i] = j 意味着原矩阵的第 i 行，应该映射到新矩阵的第 j 行。
        # 这样原 W_ica[i, j] 就会变成新矩阵的对角线元素 W_ica_star[j, j]
        W_ica_star = np.zeros_like(W_ica)
        W_ica_star[col_ind, :] = W_ica[row_ind, :]

        # --- Step 3: 归一化 (Scaling) ---
        # 严格遵守论文公式 (15): 将 W_ica_star 的每一行除以其对角线元素
        print("  Normalizing (Step 3)...")
        diag_elements = np.diag(W_ica_star).copy()

        # 防止除以极小数导致数值爆炸
        zero_mask = np.abs(diag_elements) < epsilon_val
        diag_elements[zero_mask] = epsilon_val

        # 利用 numpy 广播机制进行行归一化
        W_ica_star_normalized = W_ica_star / diag_elements[:, np.newaxis]

        # --- Step 4: 求解因果交互矩阵 B ---
        # 严格遵守论文 Step 4: B = I_N - W_ICA^*
        print("  Computing causal matrix B (Step 4)...")
        B_data = np.eye(N) - W_ica_star_normalized
        np.fill_diagonal(B_data, 0)  # 物理约束：节点不能瞬间作为自身的因果，对角线强制为0

        # --- Step 4.5: 应用物理先验修正 (Physics-Informed Prior) ---
        print("  Applying Physics-Informed Prior (Scheme 1)...")
        # 默认使用 penalty_factor=0.3，保留30%的非物理连接权重（可调节）
        B_final = apply_physics_prior(B_data, physics_adj, penalty_factor)

        # --- 结果保存 ---
        b_df = pd.DataFrame(B_final)
        b_df.to_csv(str(csv_output_path), index=False)
        print(f"  Saved B matrix to: {csv_output_path.name}")

        # --- 拓扑图可视化 ---
        G = nx.DiGraph()
        nodes = [f'{i + 1}' for i in range(N)]
        G.add_nodes_from(nodes)

        # 绘图阈值：过滤掉极弱的因果边使图像更清晰
        plot_threshold = 0.05

        for i in range(N):
            for j in range(N):
                if i == j: continue
                causal_effect = B_final[i, j]
                # 注意网络图的方向：矩阵 B_final[i, j] 表示 j 是 i 的 parent
                # 即因果方向是 j -> i
                if abs(causal_effect) > plot_threshold:
                    G.add_edge(f'{j + 1}', f'{i + 1}', weight=causal_effect)

        if not G.edges():
            print(f"  Graph has no edges (threshold={plot_threshold}). Skipping plot.")
            return

        pos = nx.spring_layout(G, k=0.8, iterations=50, seed=42)
        plt.figure(figsize=(15, 12))
        nx.draw_networkx_nodes(G, pos, node_color='#70CDBE', node_size=1000, alpha=0.5)

        edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
        max_abs_weight = max(abs(w) for w in edge_weights) if edge_weights else 1
        edge_widths = [max(0.5, abs(w) / max_abs_weight * 6) for w in edge_weights]

        nx.draw_networkx_edges(G, pos,
                               edge_color=edge_weights,
                               width=edge_widths,
                               alpha=0.7,
                               arrowstyle='-|>',
                               arrowsize=20,
                               edge_cmap=plt.cm.RdBu_r,
                               edge_vmin=-max_abs_weight,
                               edge_vmax=max_abs_weight)

        nx.draw_networkx_labels(G, pos, font_size=15, font_weight='bold', font_family='Times New Roman')

        sm = plt.cm.ScalarMappable(cmap=plt.cm.RdBu_r, norm=plt.Normalize(vmin=-max_abs_weight, vmax=max_abs_weight))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=plt.gca(), orientation='vertical', fraction=0.05, shrink=0.65, pad=0.04)
        cbar.ax.tick_params(labelsize=18)
        cbar.set_label('Causal Strength (B Matrix)', rotation=270, labelpad=20, fontsize=20)

        plt.title(f'Causal Graph for {base_name}', fontsize=18)
        plt.axis('off')

        plt.savefig(str(plot_output_path), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved graph to: {plot_output_path.name}")

    except Exception as e:
        print(f"  !! FAILED: {npz_file_path.name}. Error: {e}", file=sys.stderr)
        if 'plt' in locals():
            plt.close()

# --- 5. 主执行器 ---

def main():
    np.set_printoptions(precision=4, suppress=True)
    warnings.filterwarnings('ignore', category=UserWarning)

    print(f"--- 批量因果分析 (含物理先验) 开始 ---")

    try:
        CSV_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        PLOT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        print(f"错误: 无法创建输出目录: {e}", file=sys.stderr)
        return

    # --- (!!) NEW: 自动生成或加载物理拓扑 ---

    # 1. 检查是否存在，如果不存在则自动调用 pandapower 生成
    if not TOPOLOGY_FILE.exists():
        print("检测到物理拓扑文件不存在，正在尝试自动生成...")
        # 修改后 (传入 CASE_NAME)
        generate_edge_adjacency(TOPOLOGY_FILE, case_name=CASE_NAME)

    # 2. 预读取数据维度 (用于校验)
    npz_files = sorted(NPZ_INPUT_DIR.glob('anomaly_F*.npz'))
    if not npz_files:
        print("未找到数据文件 (anomaly_F*.npz)。")
        return

    temp_data = np.load(str(npz_files[0]))
    N_expected = temp_data['anomalyMatrix'].shape[0]
    print(f"检测到数据维度 (Features): {N_expected}")

    # 3. 加载物理拓扑
    physics_adj = load_line_topology(TOPOLOGY_FILE, N_expected)

    if physics_adj is not None:
        print("成功加载物理拓扑矩阵。方案一 (物理引导先验) 已激活。")
    else:
        # 如果维度对不上 (例如 pandapower 版本不同导致边数不同)，给予提示
        print("警告: 物理拓扑加载失败或维度不匹配。")
        print("      可能原因: 你的仿真数据中的线路数量与 pandapower 标准 IEEE39 案例不一致。")
        print("      将退化为纯数据驱动模式。")

    # --- 后续逻辑保持不变 ---
    print(f"找到 {len(npz_files)} 个 .npz 文件。开始处理...")

    for mask in MASK_VALUES:
        print(f"\n--- 处理掩码值: {mask} ---")
        mask_dir = f'mask_{mask}'
        csv_subdir = CSV_OUTPUT_DIR / mask_dir
        plot_subdir = PLOT_OUTPUT_DIR / mask_dir
        csv_subdir.mkdir(parents=True, exist_ok=True)
        plot_subdir.mkdir(parents=True, exist_ok=True)

        processed_count = 0
        for npz_file in npz_files:
            analyze_causal_matrix(npz_file, csv_subdir, plot_subdir, physics_adj=physics_adj, penalty_factor=mask)
            processed_count += 1
            if processed_count % 10 == 0:
                print(f"  (已处理 {processed_count}/{len(npz_files)})")

    print(f"\n--- 处理完毕 ---")

if __name__ == "__main__":
    main()

