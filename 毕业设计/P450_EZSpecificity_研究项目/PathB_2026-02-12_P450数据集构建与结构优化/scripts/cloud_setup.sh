#!/bin/bash
# ==============================================================================
# EZSpecificity 云服务器一键配置脚本
#
# 用法：
#   1. 从本地 scp 这个脚本到云服务器
#   2. 在云服务器上执行: bash cloud_setup.sh
#   3. 然后从本地传 .pt 数据和训练脚本
#
# 适用于: 智川云 / AutoDL 等 PyTorch 镜像的 Linux 服务器
# 要求: 已装 PyTorch + CUDA（镜像自带）
# ==============================================================================

set -e  # 出错即停

echo "=============================================="
echo "EZSpecificity Cloud Server Setup"
echo "=============================================="

# ----- 1. 检查基础环境 -----
echo ""
echo "[1/5] 检查基础环境..."
python -c "import torch; print(f'  PyTorch {torch.__version__}')"
python -c "import torch; print(f'  CUDA available: {torch.cuda.is_available()}')"
python -c "import torch; print(f'  CUDA version: {torch.version.cuda}')" 2>/dev/null || echo "  CUDA version: N/A (无卡模式)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo "  GPU: N/A (无卡模式，正常)"

# ----- 2. 安装依赖 -----
echo ""
echo "[2/5] 安装 Python 依赖..."

# 获取 PyTorch 和 CUDA 版本用于匹配 PyG 包
TORCH_VERSION=$(python -c "import torch; print(torch.__version__.split('+')[0])")
CUDA_VERSION=$(python -c "import torch; print(torch.version.cuda.replace('.', '')[:3])" 2>/dev/null || echo "124")
echo "  PyTorch: $TORCH_VERSION, CUDA: $CUDA_VERSION"

# 基础依赖
pip install -q pytorch_lightning==1.9.0 easydict rdkit-pypi lmdb tqdm matplotlib 2>&1 | tail -1
echo "  基础依赖已安装"

# PyTorch Geometric
pip install -q torch_geometric 2>&1 | tail -1
echo "  torch_geometric 已安装"

# PyG 扩展包（需要匹配 PyTorch + CUDA 版本）
# 尝试多个 wheel 源，自动匹配
PYG_WHEEL_URL="https://data.pyg.org/whl/torch-${TORCH_VERSION}+cu${CUDA_VERSION}.html"
echo "  尝试 PyG 扩展包: $PYG_WHEEL_URL"
pip install -q torch_scatter torch_sparse torch_cluster torch_spline_conv -f "$PYG_WHEEL_URL" 2>&1 | tail -1 || {
    echo "  [WARN] 精确版本匹配失败，尝试兼容版本..."
    # 回退：尝试不指定 CUDA 版本
    pip install -q torch_scatter torch_sparse torch_cluster torch_spline_conv 2>&1 | tail -1 || {
        echo "  [ERROR] PyG 扩展包安装失败，可能需要手动安装"
        echo "  参考: https://pytorch-geometric.readthedocs.io/en/latest/install/installation.html"
    }
}

# 验证安装
echo ""
echo "[3/5] 验证安装..."
python -c "
import torch
import torch_geometric
import pytorch_lightning as pl
from torch_geometric.nn import knn_graph
print(f'  torch_geometric: {torch_geometric.__version__}')
print(f'  pytorch_lightning: {pl.__version__}')
print(f'  knn_graph: OK')
print('  所有依赖验证通过!')
" || {
    echo "  [ERROR] 依赖验证失败，请检查上面的错误信息"
    exit 1
}

# ----- 3. 创建目录结构 -----
echo ""
echo "[4/5] 创建项目目录..."
PROJECT_DIR=~/EZSpecificity
mkdir -p $PROJECT_DIR/src
mkdir -p $PROJECT_DIR/scripts/10_Step10_pt训练管线
mkdir -p $PROJECT_DIR/data/10_Step10_pt训练
mkdir -p $PROJECT_DIR/results/10_Step10_pt训练/checkpoints
mkdir -p $PROJECT_DIR/logs/10_Step10_pt训练
echo "  项目目录: $PROJECT_DIR"
ls -d $PROJECT_DIR/*/

# ----- 4. 显示后续步骤 -----
echo ""
echo "[5/5] 配置完成!"
echo "=============================================="
echo ""
echo "接下来从本地执行以下命令传文件:"
echo ""
echo "# 1. 传模型源代码 (src/)"
echo "scp -r -P PORT D:\\EZSpecificity_Project\\src\\* USER@HOST:~/EZSpecificity/src/"
echo ""
echo "# 2. 传训练脚本"
echo "scp -P PORT pt_dataset.py USER@HOST:~/EZSpecificity/scripts/10_Step10_pt训练管线/"
echo "scp -P PORT main_training_pt.py USER@HOST:~/EZSpecificity/scripts/10_Step10_pt训练管线/"
echo "scp -P PORT server_config.yml USER@HOST:~/EZSpecificity/scripts/"
echo ""
echo "# 3. 传 .pt 训练数据 (57GB, 用 SFTP 客户端拖拽更方便)"
echo "# 本地路径: data/10_Step10_pt训练/ezspec_pt_v1/"
echo "# 服务器路径: ~/EZSpecificity/data/10_Step10_pt训练/ezspec_pt_v1/"
echo ""
echo "# 4. 开始训练"
echo "cd ~/EZSpecificity/src && python ../scripts/10_Step10_pt训练管线/main_training_pt.py \\"
echo "  --config ../scripts/server_config.yml \\"
echo "  --cache-dir ../data/10_Step10_pt训练/ezspec_pt_v1 \\"
echo "  --edge-mode legacy_bug --num-workers 2 --batch-size 24"
echo ""
echo "=============================================="
