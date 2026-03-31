import os
import numpy as np
from scipy.io import loadmat
from pathlib import Path
import sys
import warnings
import re

# --- 1. 配置区 ---

CASE_NAME = 'IEEE118'  # 根据你生成的系统进行修改 (如 IEEE14, IEEE39, IEEE118)
ROOT_DIR = Path(f'./data/dc/{CASE_NAME}')
OUTPUT_DIR = ROOT_DIR / 'anomaly_data_npz'

CAPACITY_VAR_NAME = 'capacityValues'
POWER_FLOW_VAR_NAME = 'powerFlowData'


def process_data():
    """
    主函数: 遍历、计算、聚合并保存异常指数矩阵。
    适配 LHS 联合采样新结构: failure_XXXX / data_LXXXX_RXXXX.mat
    """
    print(f"--- 开始计算并聚合 {CASE_NAME} 的异常指数矩阵 ---")
    print(f"--- 目标路径结构: {ROOT_DIR.name} / failure_XXXX / data_L*_R*.mat ---")

    warnings.filterwarnings('ignore', category=RuntimeWarning, message='divide by zero')
    warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value')

    try:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        print(f"错误: 无法创建输出目录: {e}", file=sys.stderr)
        return

    # 2. 加载全局容量向量
    # 尝试多种可能的命名规则，确保能抓到 linkset
    linkset_files = list((ROOT_DIR / 'linkset').glob('*.mat'))
    if not linkset_files:
        print("致命错误: 无法找到 linkset 文件夹或容量数据。", file=sys.stderr)
        return

    linkset_file_path = linkset_files[0]  # 默认取第一个
    try:
        capacity_data = loadmat(str(linkset_file_path))
        capacity_vector = capacity_data[CAPACITY_VAR_NAME].flatten()
        print(f"成功加载容量向量 (共 {capacity_vector.shape[0]} 条线路) from: {linkset_file_path.name}")
    except Exception as e:
        print(f"致命错误: 无法加载容量数据。{e}", file=sys.stderr)
        return

    # 3. 查找所有的故障文件夹
    failure_dirs = sorted(ROOT_DIR.glob('failure_*'))
    if not failure_dirs:
        print(f"警告: 在 {ROOT_DIR} 中未找到任何 'failure_*' 文件夹。", file=sys.stderr)
        return

    print(f"找到 {len(failure_dirs)} 个故障文件夹。开始处理...")
    processed_failures = 0

    # 用于解析文件名的正则表达式: 匹配 data_L0p9000_R0p5000.mat 这种格式
    pattern = re.compile(r'data_L([\dp]+)_R([\dp]+)\.mat')

    # L1 循环: 遍历【故障】 (failure_XXXX)
    for failure_dir in failure_dirs:
        if not failure_dir.is_dir():
            continue

        failure_str = failure_dir.name.split('_')[-1]

        all_anomaly_vectors = []
        load_scales = []
        renew_scales = []

        # L2 循环: 遍历该故障下的所有 LHS 样本文件
        mat_files = sorted(failure_dir.glob('data_L*_R*.mat'))

        if not mat_files:
            continue

        for filepath in mat_files:
            try:
                # 解析 L 和 R 的系数
                match = pattern.search(filepath.name)
                if match:
                    l_val = float(match.group(1).replace('p', '.'))
                    r_val = float(match.group(2).replace('p', '.'))
                else:
                    continue

                data = loadmat(str(filepath))
                power_flow_matrix = data[POWER_FLOW_VAR_NAME]

                if power_flow_matrix.shape[1] < 2:
                    continue

                # --- 核心计算: 异常指数 S_i ---
                # S_i[t] = (P_i[1] - P_i[0]) / P_i^max
                power_flow_delta_vec = power_flow_matrix[:, 1] - power_flow_matrix[:, 0]
                anomaly_vector = power_flow_delta_vec / capacity_vector

                anomaly_vector[np.isinf(anomaly_vector)] = 0
                anomaly_vector = np.nan_to_num(anomaly_vector)

                all_anomaly_vectors.append(anomaly_vector)
                load_scales.append(l_val)
                renew_scales.append(r_val)

            except Exception as e:
                print(f"  警告: 跳过 {filepath.name}: {e}", file=sys.stderr)

        # 保存聚合矩阵
        if not all_anomaly_vectors:
            continue

        # 矩阵形状: [线路数, LHS样本数]
        aggregated_matrix = np.stack(all_anomaly_vectors, axis=1)

        output_name = f"anomaly_F{failure_str}.npz"
        output_path = OUTPUT_DIR / output_name

        try:
            # 将双重变量的特征一并保存
            np.savez_compressed(
                str(output_path),
                anomalyMatrix=aggregated_matrix,
                loadScales=np.array(load_scales, dtype=np.float32),
                renewScales=np.array(renew_scales, dtype=np.float32)
            )

            processed_failures += 1
            if processed_failures % 10 == 0:
                print(f"  ({processed_failures}) 已处理: {output_name} (含 {len(load_scales)} 个双变量 LHS 样本)")

        except Exception as e:
            print(f"  !! 严重错误: 无法保存 {output_name}: {e}", file=sys.stderr)

    print(f"--- 处理完毕。共成功生成了 {processed_failures} 个聚合文件。---")


if __name__ == "__main__":
    process_data()