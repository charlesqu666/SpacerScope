# Web 入口层

该目录存放 `spacerscope-web` 独立本地 Web GUI 的应用层文件。

包含：
- `web_main.cpp`：HTTP 服务入口与路由
- `web_job.hpp/cpp`：单任务状态机、子进程管理、SSE 事件、产物下载
- `web_assets.hpp`：内嵌静态 HTML/CSS/JS 前端

说明：
- 本目录不包含搜索内核。
- `spacerscope-web` 通过子进程方式调用 `spacerscope`。
- Web v1 为单用户、单任务、无认证设计。
