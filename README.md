# ldbc_partition

## 代码目录
```
BuildScript/  建库脚本
  -buildDistributeDB.sh message分片建库
  -specialSlice.sh  forum分片建库
distribute_header/ 分布式建库特殊header
ldbcPartition.cpp 划分代码
*route/  生成路由表目录
*output/ 划分数据生成目录
  -forum/  forum分片数据路径
  -messages/  message分片数据路径

* 为执行划分后生成的数据目录
```

## 划分ldbcPartion.cpp
> N个message分片+1个forum分片，分片数目指message分片数目
### 编译
```
$ g++ ldbcPartition.cpp -o partition
```
### 运行
```
# sf0.1三分片
$ ./partition 0.1 /home/shared/data/social_network-csv_composite-longdateformatter-sf0.1/dynamic 3
# 执行后，当前目录将生成output数据文件夹(含message分片和forum分片数据),route为路由表文件夹(message2part.csv、person2messagePart.csv)
# sf30三分片
$ ./partition 30 /home/shared/data/social_network-csv_composite-longdateformatter-sf30/dynamic 3

# 参数信息
$ ./partition -h
Usage: ./program_name [SCALE_FACTOR] [DIRECTORY] [PART]

Arguments:
  SCALE_FACTOR    The ldbc scale factor value (e.g., 0.1, 3, 30, 300). Default is 0.1.
  DIRECTORY       The path to the initial ldbc dynamic directory. Default is:
                  "/home/shared/data/social_network-csv_composite-longdateformatter-sf0.1/dynamic/".
  PART            An integer specifying the message part number. Default is 3.

Options:
  -h, --help      Display this help message and exit.

Examples:
  ./program_name                         # Uses default SCALE_FACTOR, DIRECTORY, and PART
  ./program_name 0.1                       # Uses SCALE_FACTOR = 1, default DIRECTORY and PART
  ./program_name 0.1 /path/to/dir          # Uses SCALE_FACTOR = 1, DIRECTORY = /path/to/dir, and default PART
  ./program_name 0.1 /path/to/dir 5        # Uses SCALE_FACTOR = 1, DIRECTORY = /path/to/dir, and PART = 5
```
## 建库 BuildScript
```
# 在gpstore目录执行
# forum分片
$ sh BuildScript/specialSlice.sh
# message分片，输入分片id参数
$ sh BuildScript/buildDistributeDB.sh id
```
### 脚本配置
```
# ldbc原始数据目录
ldbcPath="/home/huzheyuan/ldbc_partition/input/social_network-csv_composite-longdateformatter-sf0.1"
# gpstore bin目录
binPath="/home/huzheyuan/latest-build/gpstore/bin"
# 分布式特殊header目录
divideHeader="/home/huzheyuan/ldbc_partition/distribute_header"
# 分布式数据(forum分片/message分片数据目录)
DataRoot="/home/huzheyuan/ldbc_partition/output/forum"
# 建库目录
DataBaseHome="./SF0.1_Divide3_3_end"
```