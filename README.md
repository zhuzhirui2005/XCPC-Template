# XCPC-Template

XCPC 算法模板，以及自动生成的 [C++ 头文件 PDF](headers.pdf)。

## OI-Wiki PDF

[下载本仓库自动同步的最新版 OI-wiki.pdf](https://github.com/zhuzhirui2005/XCPC-Template/releases/download/oi-wiki-latest/OI-wiki.pdf)

`Update OI-Wiki PDF` 每天 UTC 01:00 检查 [OI-Wiki-export 的最新 Release](https://github.com/OI-wiki/OI-Wiki-export/releases/latest)，下载并校验 PDF，然后更新本仓库 `oi-wiki-latest` Release 的附件。也可以在 Actions 页面手动运行。

PDF 已超过 GitHub 普通 Git 文件的 100 MiB 限制，因此新版改用 Release 附件分发，不再提交二进制文件，也不需要 Git LFS。

仓库根目录现有的 `OI-wiki.pdf` 保留为历史快照，不再自动更新；请使用上方链接获取新版。首次发布完成前，也可从上述上游 Release 下载。下载后可以将 PDF 保存在本地以便离线使用。

修改此 workflow 的 Pull Request 会执行下载和完整性校验，但不会发布 Release。合并到 `master` 后会自动执行首次发布。
