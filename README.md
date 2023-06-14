# AQP

## 简介

用 C 实现的 aqp，可以在规定时间内求得准确值（卷翻天）。

## 运行须知

- 仅支持 Linux 系统（测试环境为 Ubuntu 20.04 和 dlab）。
- 请保证 `aqp.py` 和 `codes/` 的相对路径关系，`aqp.py` 的运行依赖于 `2021201626/codes/` 下的 `kdtree_aqp.py` 和 `libaqp.so`。
- 程序运行需要足够大 (上传的版本需要4GB) 的临时硬盘空间用作数据目录，默认是 `/2021201626/tmp` 。**如果不满足这个条件**，请在 `kdtree_aqp.py` 中手动修改源码中的 `DATA_DIR` 常量为一个可用的绝对路径（支持自动创建目录）。

## 手动编译（如果要在其他平台上运行）

Linux 下直接在 `2021201626/codes/` 下执行 `make` 即可编译出 `libaqp.so`。
如果要在其他平台上运行请参考 `Makefile` 文件。同时在 `kdtree_aqp.py` 中查找字段 `lib = CDLL(osp.join(CODE_DIR,'libaqp.so'))`，将其中 `osp.join(CODE_DIR,'libaqp.so')` 改为对应的动态库所在路径。

> 通常情况下，无论什么平台，你只需要重新编译一次 `g++ -Ofast -shared -fPIC -o libaqp.so libaqp.cc` 即可

## 报告

实验报告见 [docs/report.pdf](./docs/report.pdf)
