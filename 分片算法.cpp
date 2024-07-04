struct Vertex {
    int cnt;    //一跳邻居个数
    set 2-hop;  //二跳内的邻居集合
    int seed = -1;   //入选分片
};

int partition[];  //每个点的最新分片状态

int vertex_cnt[];   //每个分片的点个数
int threshold;    //每个分片的最大点个数
int p = 0;        //分片个数


Code:

while(存在点K，使得(K.seed == -1))  //存在点未入选种子集
{
    取点K， 使得(K.seed == -1 && K.cnt最大)
    
    
    //将K及其两跳邻居加入p分片
    K.seed = ++ p;
    partition[K] = p;
    vertex_cnt[p] ++;
    
    for(U : K.2-hop)
    {
        partition[U] = p;
        vertex_cnt[p] ++;
    }
    
    
    //尝试将K的一跳邻居L加入p分片的种子集中
    while(vertex_cnt[p] < threshold)
    {
        取点L, 使得(L.seed == -1 && L与K相邻 && L.cnt最小)
        L.seed = p;
        
        //将L的两跳邻居加入p分片
        for(U : L.2-hop)
            if(partition[U] < p)    //该点不在p分片里
            {
                partition[U] = p;
                vertex_cnt[p] ++;
            }
    }
    //该算法将导致每个分片都会超出threshold一点，但可以避免多次重复搜索同一个点的两跳邻居
    //可以考虑将threshold适当缩小
    //也可以选择先计算L.2-hop与p分片中不重叠的点个数，加上vertex_cnt[p]再与threshold做比较，缺点是可能会导致多次计算同一个节点的2-hop是否重叠
}
