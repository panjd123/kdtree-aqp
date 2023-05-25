# AQP 2021201626

重要！！

- 仅支持 Linux 系统（测试环境为 Ubuntu 20.04 和 dlab）。
- 请保证 `aqp.py` 和 `codes/` 的相对路径关系，`aqp.py` 的运行依赖于 `2021201626/codes/` 下的 `kdtree_aqp.py` 和 `libaqp.so`。
- 程序运行需要足够大 (上传的版本需要4GB) 的临时硬盘空间用作数据目录，默认是 `/home/mw/temp` 。**如果不满足这个条件**，请在 `kdtree_aqp.py` 中手动修改源码中的 `DATA_DIR` 常量为一个可用的绝对路径（支持自动创建目录）。

> 根据 `dlab` 的文档，`/home/mw/temp` 是指定的临时数据存放区，容量应当支持至少 100G，但如果同时测试且不清除之前的人的数据可能会容量不足。

## 手动编译（如果要在其他平台上运行）

Linux 下直接在 `2021201626/codes/` 下执行 `make` 即可编译出 `libaqp.so`。
如果要在其他平台上运行请参考 `Makefile` 文件。同时在 `kdtree_aqp.py` 中
查找字段 `lib = CDLL(osp.join(CODE_DIR,'libaqp.so'))`，将其中 `osp.join(CODE_DIR,'libaqp.so')` 改为对应的动态库所在路径。