# GitHub 发布指南

本文档指导您如何将项目发布到 GitHub 并进行管理。

## 📋 发布前检查清单

在发布到 GitHub 之前，确保完成以下步骤：

- [x] ✅ 创建 `.gitignore` 文件
- [x] ✅ 编写完整的 `README.md`
- [x] ✅ 添加 `LICENSE` 文件
- [x] ✅ 创建 `requirements.txt`
- [x] ✅ 编写 `CHANGELOG.md`
- [x] ✅ 添加 `CONTRIBUTING.md`
- [ ] ⚠️ 清理临时文件和敏感数据
- [ ] ⚠️ 测试所有功能
- [ ] ⚠️ 更新 README 中的 GitHub 链接

## 🚀 发布步骤

### 步骤 1：初始化 Git 仓库

```bash
# 进入项目目录
cd "F:\2025-08-赖草根系微生物多组学项目\2025-11-培养基微生物16S重新检验结果"

# 初始化 Git 仓库
git init

# 添加所有文件
git add .

# 查看状态（确认要提交的文件）
git status

# 首次提交
git commit -m "Initial commit: Sanger 16S rRNA Analysis Toolkit v1.0.0"
```

### 步骤 2：在 GitHub 上创建仓库

1. **登录 GitHub**
   - 访问 https://github.com
   - 登录您的账号

2. **创建新仓库**
   - 点击右上角的 "+" → "New repository"
   - 填写仓库信息：
     - **Repository name**: `sanger-16s-analysis`（或您喜欢的名称）
     - **Description**: `智能的 Sanger 16S rRNA 测序数据分析工具包`
     - **Visibility**: Public（公开）或 Private（私有）
     - **不要**勾选 "Initialize this repository with a README"（我们已有 README）

3. **创建仓库**
   - 点击 "Create repository"

### 步骤 3：连接本地仓库到 GitHub

```bash
# 添加远程仓库（替换 YOUR_USERNAME 为您的 GitHub 用户名）
git remote add origin https://github.com/YOUR_USERNAME/sanger-16s-analysis.git

# 验证远程仓库
git remote -v

# 推送到 GitHub（首次推送）
git branch -M main
git push -u origin main
```

### 步骤 4：更新 README 中的链接

在 `README.md` 中，将所有 `YOUR_USERNAME` 替换为您的实际 GitHub 用户名：

```bash
# 使用文本编辑器批量替换
# 或使用命令行（Linux/macOS）:
sed -i 's/YOUR_USERNAME/your_actual_username/g' README.md
sed -i 's/YOUR_USERNAME/your_actual_username/g' CONTRIBUTING.md
sed -i 's/YOUR_USERNAME/your_actual_username/g' CHANGELOG.md
sed -i 's/YOUR_USERNAME/your_actual_username/g' examples/README.md
sed -i 's/YOUR_USERNAME/your_actual_username/g' docs/PROJECT_STRUCTURE.md

# Windows PowerShell:
(Get-Content README.md) -replace 'YOUR_USERNAME', 'your_actual_username' | Set-Content README.md
(Get-Content CONTRIBUTING.md) -replace 'YOUR_USERNAME', 'your_actual_username' | Set-Content CONTRIBUTING.md
(Get-Content CHANGELOG.md) -replace 'YOUR_USERNAME', 'your_actual_username' | Set-Content CHANGELOG.md
(Get-Content examples/README.md) -replace 'YOUR_USERNAME', 'your_actual_username' | Set-Content examples/README.md
(Get-Content docs/PROJECT_STRUCTURE.md) -replace 'YOUR_USERNAME', 'your_actual_username' | Set-Content docs/PROJECT_STRUCTURE.md
```

然后提交更改：

```bash
git add .
git commit -m "docs: 更新 GitHub 用户名链接"
git push
```

### 步骤 5：添加仓库主题和描述

在 GitHub 仓库页面：

1. **添加主题标签（Topics）**
   - 点击仓库名称下方的 "Add topics"
   - 添加相关标签：
     - `bioinformatics`
     - `16s-rrna`
     - `sanger-sequencing`
     - `sequence-analysis`
     - `blast`
     - `python`
     - `genomics`
     - `microbiology`

2. **编辑仓库描述**
   - 点击右侧的 ⚙️ 图标
   - 添加简短描述和网站链接（如果有）

### 步骤 6：配置仓库设置

在仓库的 Settings 中：

1. **Features**
   - ✅ Issues（启用问题追踪）
   - ✅ Discussions（启用讨论功能）
   - ✅ Projects（可选）
   - ✅ Wiki（可选）

2. **Pull Requests**
   - ✅ Allow merge commits
   - ✅ Allow squash merging
   - ✅ Allow rebase merging

3. **Branches**
   - 设置 `main` 为默认分支
   - 可选：添加分支保护规则

## 📝 创建 Release

### 创建第一个版本发布

1. **在 GitHub 上创建 Release**
   - 进入仓库页面
   - 点击右侧的 "Releases" → "Create a new release"

2. **填写 Release 信息**
   - **Tag version**: `v1.0.0`
   - **Release title**: `v1.0.0 - 初始发布`
   - **Description**: 从 CHANGELOG.md 复制内容
   
   ```markdown
   ## 🎉 首次发布
   
   ### 主要特性
   - ✨ 智能文件名识别系统
   - ✨ 完整的分析流程
   - ✨ BLAST 断点续传
   - ✨ 交互式和批处理模式
   - ✨ 灵活的配置系统
   
   ### 核心模块
   - `filename_parser.py` - 文件名解析
   - `graph_based_assembly.py` - 序列组装
   - `ncbi_blast_analysis.py` - BLAST 分析
   - `sanger_16s_pipeline.py` - 流程管理
   
   ### 安装
   ```bash
   pip install -r requirements.txt
   ```
   
   ### 快速开始
   ```bash
   python sanger_16s_pipeline.py
   ```
   
   详细文档请查看 [README.md](README.md)
   ```

3. **发布**
   - 点击 "Publish release"

## 🎨 美化仓库

### 添加徽章（Badges）

README 中已包含基础徽章，您还可以添加：

```markdown
<!-- 代码质量 -->
![Code Quality](https://img.shields.io/badge/code%20quality-A-brightgreen)

<!-- 文档 -->
![Documentation](https://img.shields.io/badge/docs-passing-brightgreen)

<!-- 活跃度 -->
![Maintenance](https://img.shields.io/badge/maintained-yes-brightgreen)

<!-- 下载量（需要实际数据）-->
![Downloads](https://img.shields.io/github/downloads/YOUR_USERNAME/sanger-16s-analysis/total)

<!-- Stars -->
![GitHub stars](https://img.shields.io/github/stars/YOUR_USERNAME/sanger-16s-analysis?style=social)
```

### 添加 GitHub Actions（可选）

创建 `.github/workflows/tests.yml`：

```yaml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: [3.7, 3.8, 3.9, '3.10']
    
    steps:
    - uses: actions/checkout@v2
    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v2
      with:
        python-version: ${{ matrix.python-version }}
    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install -r requirements.txt
    - name: Run tests
      run: |
        python -m pytest tests/
```

## 📢 推广项目

### 1. 完善项目描述

确保 README 包含：
- 清晰的项目介绍
- 吸引人的特性列表
- 简单的快速开始指南
- 详细的使用文档
- 贡献指南

### 2. 添加截图和 GIF

如果有 GUI 或命令行界面，添加截图：

```markdown
## 📸 截图

### 交互式模式
![Interactive Mode](docs/images/interactive_mode.png)

### 分析结果
![Results](docs/images/results.png)
```

### 3. 分享到社区

- 在相关的生物信息学论坛分享
- 在 Twitter/微博上发布
- 在学术社交网络（ResearchGate）分享
- 在相关的 Reddit 子版块发布

### 4. 撰写博客文章

写一篇详细的介绍文章：
- 项目背景和动机
- 主要功能和特性
- 使用教程
- 技术实现细节

## 🔄 日常维护

### 提交代码规范

```bash
# 功能开发
git checkout -b feature/new-feature
# ... 进行修改 ...
git add .
git commit -m "feat: 添加新功能描述"
git push origin feature/new-feature
# 在 GitHub 上创建 Pull Request

# Bug 修复
git checkout -b fix/bug-description
# ... 修复 bug ...
git add .
git commit -m "fix: 修复 bug 描述"
git push origin fix/bug-description

# 文档更新
git checkout -b docs/update-readme
# ... 更新文档 ...
git add .
git commit -m "docs: 更新文档"
git push origin docs/update-readme
```

### 处理 Issues

1. **及时响应**
   - 尽快回复新的 Issue
   - 感谢用户的反馈

2. **标签管理**
   - `bug` - Bug 报告
   - `enhancement` - 功能建议
   - `documentation` - 文档相关
   - `good first issue` - 适合新手
   - `help wanted` - 需要帮助

3. **关闭 Issue**
   - 修复后关联提交：`Fixes #123`
   - 说明解决方案

### 审查 Pull Requests

1. **代码审查**
   - 检查代码质量
   - 确保符合规范
   - 测试功能

2. **反馈**
   - 提供建设性意见
   - 感谢贡献者

3. **合并**
   - 使用合适的合并策略
   - 更新 CHANGELOG

## 📊 项目统计

### 使用 GitHub Insights

查看项目统计：
- **Traffic**: 访问量和克隆数
- **Contributors**: 贡献者列表
- **Community**: 社区健康度
- **Pulse**: 项目活跃度

### 添加统计徽章

```markdown
![GitHub contributors](https://img.shields.io/github/contributors/YOUR_USERNAME/sanger-16s-analysis)
![GitHub last commit](https://img.shields.io/github/last-commit/YOUR_USERNAME/sanger-16s-analysis)
![GitHub issues](https://img.shields.io/github/issues/YOUR_USERNAME/sanger-16s-analysis)
![GitHub pull requests](https://img.shields.io/github/issues-pr/YOUR_USERNAME/sanger-16s-analysis)
```

## 🎯 最佳实践

### 1. 保持代码整洁
- 定期重构
- 遵循编码规范
- 添加注释和文档

### 2. 及时更新文档
- 功能变更时更新 README
- 维护 CHANGELOG
- 更新示例代码

### 3. 版本管理
- 遵循语义化版本
- 创建 Git 标签
- 发布 Release

### 4. 社区互动
- 回复 Issues
- 审查 Pull Requests
- 感谢贡献者

### 5. 持续改进
- 收集用户反馈
- 修复 Bug
- 添加新功能

## 🔒 安全注意事项

### 不要提交的内容

- ❌ 个人敏感信息
- ❌ API 密钥和密码
- ❌ 大型数据文件（使用 Git LFS）
- ❌ 临时文件和缓存
- ❌ 编译产物

### 使用 .gitignore

确保 `.gitignore` 包含：
```
# Python
__pycache__/
*.py[cod]
*.egg-info/
venv/

# 数据文件
*.ab1
*.fasta
*.seq
*.pdf

# 结果文件
*_results/
blast_results/

# 日志
*.log

# 临时文件
temp/
*.tmp
```

## 📚 相关资源

- [GitHub 文档](https://docs.github.com/)
- [Git 教程](https://git-scm.com/book/zh/v2)
- [语义化版本](https://semver.org/lang/zh-CN/)
- [Keep a Changelog](https://keepachangelog.com/zh-CN/)
- [开源许可证选择](https://choosealicense.com/)

## ❓ 常见问题

### Q: 如何更改仓库名称？
A: Settings → Repository name → Rename

### Q: 如何删除敏感信息的提交历史？
A: 使用 `git filter-branch` 或 BFG Repo-Cleaner

### Q: 如何处理大文件？
A: 使用 Git LFS（Large File Storage）

### Q: 如何设置协作者？
A: Settings → Collaborators → Add people

### Q: 如何创建组织仓库？
A: 在组织页面创建仓库，而不是个人账户

---

祝您的项目在 GitHub 上获得成功！🎉

如有问题，请参考 [GitHub 官方文档](https://docs.github.com/) 或提交 Issue。

