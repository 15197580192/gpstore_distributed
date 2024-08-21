# sh BuildGraph/specialSlice.sh 

# ldbc原始数据
ldbcPath="/home/huzheyuan/ldbc_partition/input/social_network-csv_composite-longdateformatter-sf0.1"
dynamicData="$ldbcPath/dynamic"
staticData="$ldbcPath/static"
dynamicHeader="$ldbcPath/header/dynamic"
staticHeader="$ldbcPath/header/static"
# gpstore bin目录
binPath="/home/huzheyuan/latest-build/gpstore/bin"
# 分布式特殊header
divideHeader="/home/huzheyuan/ldbc_partition/distribute_header"
# 分布式数据(forum分片)
DataRoot="/home/huzheyuan/ldbc_partition/output/forum"
# 建库目录
DataBaseHome="./SF0.1_Divide3_3_end"

$binPath/gdrop -db "$DataBaseHome/db_forum_container_of"
$binPath/gpbuild -db "$DataBaseHome/db_forum_container_of"\
    --nodes=Forum=$dynamicHeader/Forum.csv,$dynamicData/forum_0_0.csv\
    --nodes=Message,Post=$divideHeader/Post-forumPart.csv,$DataRoot/post_0_0.csv\
    --nodes=Person=$dynamicHeader/Person.csv,$dynamicData/person_0_0.csv\
    --relationships=CONTAINER_OF=$dynamicHeader/Forum_containerOf_Post.csv,$dynamicData/forum_containerOf_post_0_0.csv\
    --relationships=HAS_MODERATOR=$dynamicHeader/Forum_hasModerator_Person.csv,$dynamicData/forum_hasModerator_person_0_0.csv\
    --relationships=HAS_MEMBER=$dynamicHeader/Forum_hasMember_Person.csv,$dynamicData/forum_hasMember_person_0_0.csv\
    --relationships=POST_HAS_CREATOR=$dynamicHeader/Post_hasCreator_Person.csv,$dynamicData/post_hasCreator_person_0_0.csv