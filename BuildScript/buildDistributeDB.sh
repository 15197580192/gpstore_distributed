#!/bin/bash
#输入参数是分片id
#sh BuildScript/buildDistributeDB.sh

# ldbc原始数据
ldbcPath="/home/huzheyuan/ldbc_partition/input/social_network-csv_composite-longdateformatter-sf0.1"
HeaderPath="$ldbcPath/headers"
dynamicData="$ldbcPath/dynamic"
dynamicHeader="$ldbcPath/header/dynamic"
staticHeader="$ldbcPath/header/static"
# 分布式特殊header
divideHeader="/home/huzheyuan/ldbc_partition/distribute_header"
# gpstore bin目录
binPath="/home/huzheyuan/latest-build/gpstore/bin"
# 分布式数据(message分片)
DataRoot="/home/huzheyuan/ldbc_partition/output/messages"
DataPath="$DataRoot/$1"
# 建库目录
DataBaseHome="./SF0.1_Divide3_3_end"
databaseName="$DataBaseHome/db_$1"

$binPath/gdrop -db $databaseName
$binPath/gpbuild -db $databaseName\
        --nodes=Message,Comment=$dynamicHeader/Comment.csv,$DataPath/comment_0_0.csv\
        --nodes=Message,Comment=$divideHeader/Comment.csv,$DataPath/comment_0_0_added.csv\
        --nodes=Forum=$divideHeader/Forum.csv,$DataPath/forum_0_0.csv\
        --nodes=Message,Post=$dynamicHeader/Post.csv,$DataPath/post_0_0.csv\
        --nodes=Message,Post=$divideHeader/Post.csv,$DataPath/post_0_0_added.csv\
        --nodes=Organisation=$staticHeader/Organisation.csv,$staticData/organisation_0_0.csv\
        --nodes=Person=$dynamicHeader/Person.csv,$dynamicData/person_0_0.csv\
        --nodes=TagClass=$staticHeader/TagClass.csv,$staticData/tagclass_0_0.csv\
        --nodes=Tag=$staticHeader/Tag.csv,$staticData/tag_0_0.csv\
        --nodes=Place=$staticHeader/Place.csv,$staticData/place_0_0.csv\
        --relationships=COMMENT_HAS_CREATOR=$dynamicHeader/Comment_hasCreator_Person.csv,$DataPath/comment_hasCreator_person_0_0.csv\
        --relationships=POST_HAS_CREATOR=$dynamicHeader/Post_hasCreator_Person.csv,$DataPath/post_hasCreator_person_0_0.csv\
        --relationships=COMMENT_IS_LOCATED_IN=$dynamicHeader/Comment_isLocatedIn_Country.csv,$DataPath/comment_isLocatedIn_place_0_0.csv\
        --relationships=REPLY_OF_COMMENT=$dynamicHeader/Comment_replyOf_Comment.csv,$DataPath/comment_replyOf_comment_0_0.csv\
        --relationships=REPLY_OF_POST=$dynamicHeader/Comment_replyOf_Post.csv,$DataPath/comment_replyOf_post_0_0.csv\
        --relationships=REPLY_OF_POST_END=$dynamicHeader/Comment_replyOf_Post.csv,$DataPath/comment_replyOf_post_0_0_end.csv\
        --relationships=CONTAINER_OF=$dynamicHeader/Forum_containerOf_Post.csv,$DataPath/forum_containerOf_post_0_0.csv\
        --relationships=COMMENT_HAS_TAG=$dynamicHeader/Comment_hasTag_Tag.csv,$DataPath/comment_hasTag_tag_0_0.csv\
        --relationships=POST_HAS_TAG=$dynamicHeader/Post_hasTag_Tag.csv,$DataPath/post_hasTag_tag_0_0.csv\
        --relationships=POST_IS_LOCATED_IN=$dynamicHeader/Post_isLocatedIn_Country.csv,$DataPath/post_isLocatedIn_place_0_0.csv\
        --relationships=LIKES=$dynamicHeader/Person_likes_Comment.csv,$DataPath/person_likes_comment_0_0.csv\
        --relationships=LIKES=$dynamicHeader/Person_likes_Post.csv,$DataPath/person_likes_post_0_0.csv\
        --relationships=ORGANISATION_IS_LOCATED_IN=$staticHeader/Organisation_isLocatedIn_Place.csv,$staticData/organisation_isLocatedIn_place_0_0.csv\
        --relationships=PERSON_IS_LOCATED_IN=$dynamicHeader/Person_isLocatedIn_City.csv,$dynamicData/person_isLocatedIn_place_0_0.csv\
        --relationships=KNOWS=$dynamicHeader/Person_knows_Person.csv,$dynamicData/person_knows_person_0_0.csv\
        --relationships=STUDY_AT=$dynamicHeader/Person_studyAt_University.csv,$dynamicData/person_studyAt_organisation_0_0.csv\
        --relationships=WORK_AT=$dynamicHeader/Person_workAt_Company.csv,$dynamicData/person_workAt_organisation_0_0.csv\
        --relationships=IS_PART_OF=$staticHeader/Place_isPartOf_Place.csv,$staticData/place_isPartOf_place_0_0.csv\
        --relationships=HAS_TYPE=$staticHeader/Tag_hasType_TagClass.csv,$staticData/tag_hasType_tagclass_0_0.csv\
        --relationships=IS_SUBCLASS_OF=$staticHeader/TagClass_isSubclassOf_TagClass.csv,$staticData/tagclass_isSubclassOf_tagclass_0_0.csv\
        --relationships=HAS_INTEREST=$dynamicHeader/Person_hasInterest_Tag.csv,$dynamicData/person_hasInterest_tag_0_0.csv