# XCPC-Template

XCPC 算法模板，以及自动生成的 [C++ 头文件 PDF](headers.pdf)。

## OI-Wiki PDF

[下载本仓库自动同步的最新版 OI-wiki.pdf](https://github.com/zhuzhirui2005/XCPC-Template/releases/download/oi-wiki-latest/OI-wiki.pdf) · [上游 Release](https://github.com/OI-wiki/OI-Wiki-export/releases/latest)

`Update OI-Wiki PDF` 每天 UTC 01:23 检查上游。先比较文件大小和 SHA-256，内容未变时不下载 PDF；有更新时下载、校验，再更新 `oi-wiki-latest` Release 附件。手动运行和修改同步流程的 PR 始终执行实际下载校验。

新版 PDF 超过普通 Git 单文件的 100 MiB 限制，因此使用 Release 附件分发，不引入 Git LFS。根目录旧 `OI-wiki.pdf` 已移除；可以下载到本地离线使用，`.gitignore` 会避免误提交。删除工作区文件不清除旧 Git 历史，本次没有重写历史或强制推送。

## 自动化维护

| 工作流 | 触发方式 | 检查或产物 |
| --- | --- | --- |
| C++20 and Workflow Checks | 头文件/CI 文件的 push、PR；每周一 UTC 03:17；手动 | GCC 13/14 语法检查、actionlint、脚本回归测试 |
| Generate C++ Headers PDF | 头文件或 PDF 构建/发布脚本变更的 push、PR；手动 | 生成 PDF，Actions 中保留 7 天预览产物；仅默认分支更新根目录 `headers.pdf` |
| Update OI-Wiki PDF | 每天 UTC 01:23；同步脚本变更的 push、PR；手动 | 校验并同步 Release PDF；用于跨 job 发布的临时产物保留 1 天 |

定时触发时间是 UTC 计划时间，不保证准点开始。编译流程的文件名仍为 `daily-compile.yml`，保留原工作流历史，但不再每天重复检查相同源码。

所有流程默认只读、固定 Ubuntu 24.04、配置超时和并发控制。构建 job 不持有仓库写权限，只有默认分支的独立发布 job 可以回写/发布；PR 和非默认分支不会发布。Action 使用完整 commit SHA 固定，由 Dependabot 每月合并成一组更新 PR，不自动合并。

头文件属于需要公共前置声明和依赖的竞赛代码片段。编译脚本只在临时副本中补充这些依赖，不修改 `.hpp` 原文件；会汇总所有头文件的错误，不在第一个失败处停止。当前是 C++20 **语法检查**，不代表模板实例化后的功能、算法正确性或复杂度测试。

PDF 不再每日重装 LaTeX 重建。构建规范化换行、编码、PDF 时间戳和文档 ID，不修改源码；工具链版本变化仍可能导致二进制差异。发布前检查构建输入是否仍与分支一致，使用预期分支 SHA 做原子提交；旧构建不会覆盖新源码的 PDF，相同 PDF 不产生提交，也不使用 force push。若分支保护禁止机器人直接提交，发布会明确失败，不会绕过保护；构建产物仍可下载。

## 本地运行

需要 Python 3.10+；语法检查需要支持 C++20 的 GCC；脚本回归测试还需要 Node.js 22+。

```sh
python3 .github/scripts/headers.py check --compiler g++ --jobs 2
python3 -m unittest discover -s .github/tests -p 'test_*.py' -v
node --test .github/tests/*.cjs
```

生成 PDF 需要 Pandoc、XeLaTeX、`xeCJK`、`fvextra`、Noto CJK 和 DejaVu 字体。Ubuntu 安装命令见 PDF workflow。

```sh
python3 .github/scripts/headers.py build --output build/headers.pdf
go run github.com/rhysd/actionlint/cmd/actionlint@v1.7.12 -color
```
