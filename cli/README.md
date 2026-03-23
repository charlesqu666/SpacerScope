# CLI 入口层

该目录存放 `spacerscope` 命令行主程序的应用层文件。

包含：
- `main.cpp`：CLI 模式解析、完整流程组装
- `progress.hpp/cpp`：文本/JSONL 进度输出协议
- `postprocess.hpp/cpp`：结果评分与后处理
- `report_html.hpp/cpp`：HTML 报告生成
- `report_logo.hpp`：报告内嵌 logo

说明：
- 搜索内核不在本目录，而在上一级共享核心源码中。
- `spacerscope-web` 也会调用 `spacerscope`，但不直接链接这些应用层文件。
