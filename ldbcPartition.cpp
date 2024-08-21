#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <set>
#include <sstream>
#include <chrono>
#include <algorithm>
#include <iterator>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <sys/resource.h>
#include <unistd.h>

#include <sys/stat.h>

using namespace std;

void GetCurrentProcessMemoryUsage() {
    std::ifstream statusFile("/proc/self/status");
    std::string line;

    if (!statusFile.is_open()) {
        std::cerr << "Could not open /proc/self/status" << std::endl;
        return;
    }

    while (std::getline(statusFile, line)) {
        if (line.find("VmSize:") == 0 || line.find("VmRSS:") == 0) {
            // 提取内存大小并转换为 MB
            size_t colonPos = line.find(':');
            size_t sizeStart = colonPos + 1;
            size_t sizeEnd = line.find('k', sizeStart);
            std::string sizeStr = line.substr(sizeStart, sizeEnd - sizeStart);
            size_t memorySizeKB = std::stoul(sizeStr);
            double memorySizeMB = memorySizeKB / 1024.0; // 转换为 MB

            // 输出结果
            if (line.find("VmSize:") == 0) {
                std::cout << "Virtual Memory Size: " << memorySizeMB << " MB" << std::endl;
            } else if (line.find("VmRSS:") == 0) {
                std::cout << "Resident Set Size: (Actual memory)" << memorySizeMB << " MB" << std::endl;
            }
        }
    }
}


// 用于社区发现
struct Edge {
    string from;
    string to;
    double weight;
};

struct Point {
    double value; // 权值
    std::string name; // 其他属性（可以根据需要添加）

    // 重载小于运算符，用于堆的比较
    bool operator<(const Point& other) const {
        return value < other.value; // 小于运算符，用于最大堆
    }
};
// 定义点权值最大堆
class PointHeap {
public:
    void push(const Point& point) {
        maxHeap.push(point);
    }

    Point pop() {
        Point topPoint = maxHeap.top();
        maxHeap.pop();
        return topPoint;
    }

    Point peek() const {
        return maxHeap.top();
    }

    bool isEmpty() const {
        return maxHeap.empty();
    }

    size_t size() const {
        return maxHeap.size(); // 返回堆的大小
    }

private:
    std::priority_queue<Point> maxHeap; // 最大堆
};

class Graph {
public:
    Graph() : totalWeight(0.0) {}

    void addEdge(const string& u, const string& v, double weight) {
        bool edgeFound = false;
        for (auto& edge : adjList[u]) {
            if (edge.to == v) {
                edge.weight += weight;
                edgeFound = true;
                break;
            }
        }
        if (!edgeFound) {
            adjList[u].push_back({u, v, weight});
        }

        edgeFound = false;
        for (auto& edge : adjList[v]) {
            if (edge.to == u) {
                edge.weight += weight;
                edgeFound = true;
                break;
            }
        }
        if (!edgeFound) {
            adjList[v].push_back({v, u, weight});
        }

        totalWeight += weight;
    }

    const vector<Edge>& getEdges(const string& node) const {
        static const vector<Edge> empty;
        auto it = adjList.find(node);
        return it == adjList.end() ? empty : it->second;
    }
    vector<string> getNeighbor(const string& node) const {
        vector<string> empty;
        for (const auto& edge : getEdges(node)) {
            empty.push_back(edge.to);
        }
        return empty;
    }
    const string getMaxNeighbor(const string& node) const {
        string maxNode;
        double maxWeight = 0.0;
        for(auto neighbor: getNeighbor(node)) {
            if(getNodeWeight(neighbor) > maxWeight) {
                maxWeight = getNodeWeight(neighbor);
                maxNode = neighbor;
            }
        }
        return maxNode; // 返回的是最大邻居或者自己
    }
    const string& getFirstNeighbor(const string& node) const {
        int i=0;
        while(getNeighbor(node)[i++]!=node) return getNeighbor(node)[i];
    }

    const unordered_map<string, vector<Edge>>& getAdjList() const {
        return adjList;
    }

    double getTotalWeight() const {
        return totalWeight;
    }
    void removeNodes(const unordered_set<string>& nodesToRemove) {
        for (const auto& node : nodesToRemove) {
            totalWeight -= getNodeWeight(node);
            adjList.erase(node);
        }
        for (auto& kv : adjList) {
            kv.second.erase(remove_if(kv.second.begin(), kv.second.end(),
                                      [&nodesToRemove](const Edge& e) {
                                          return nodesToRemove.count(e.to) > 0;
                                      }),
                            kv.second.end());
        }
    }
    void print() {
        for (const auto& kv : adjList) {
            cout << kv.first << ": ";
            for (const Edge& e : kv.second) {
                cout << e.to << "(" << e.weight << ") ";
            }
            cout << endl;
        }
    }
    string getMaxNode() {
        string maxNode;
        double maxWeight = 0.0;
        for (const auto& kv : adjList) {
            double weight = getNodeWeight(kv.first);
            if (weight > maxWeight) {
                maxWeight = weight;
                maxNode = kv.first;
            }
        }
        return maxNode;
    }
    double getNodeWeight(const string& node) const {
        double weight = 0.0;
        for (const auto& edge : getEdges(node)) {
            weight += edge.weight;
        }
        return weight;
    }

private:
    unordered_map<string, vector<Edge>> adjList;
    double totalWeight;
};


// 创建文件夹路径
bool createDirectory(const std::string& path) {
    // 使用 mkdir 函数创建文件夹路径
    int status = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (status == 0) {
        return true; // 创建成功
    } else {
        return false; // 创建失败
    }
}

void dfs(std::string current, std::unordered_map<std::string, std::unordered_set<std::string>> &reply, std::vector<std::string> &tmp, std::unordered_map<std::string,std::vector<std::string>> &commentPath) {
    // 需要 from neighbor end post
    // commentPath----(fromComment neighborComment endComment)
    if(tmp.size()<=2) {
        // 最多存三个，用于合并后续邻居(起点，一步邻居，末端comment)
        tmp.push_back(current);
    }
    if(reply[current].size()==0) {
        if(tmp.size()==1) {
            // 顶点没有邻居
        } else if(tmp.size()==2) {
            // 顶点只有一步邻居
        } else {
            // 顶点有多步邻居
            tmp.pop_back();
            tmp.push_back(current);
        }
        commentPath[tmp[0]] = tmp;
        return;
    }
    for (auto i : reply[current]) {
        dfs(i, reply, tmp, commentPath);
    }
}


// dfs得到所有的message链
void dfs_c(std::string current, std::unordered_map<std::string, std::unordered_set<std::string>> &reply1, std::set<std::string> &tmp, std::unordered_map<std::string,std::unordered_map<std::string,int>> visited) {
    for (auto i : reply1[current]) {
        if(visited[current][i]==1) {
            continue;
        }
        visited[current][i] = 1;
        tmp.insert(i);
        dfs_c(i, reply1, tmp,visited);
    }
}

void loadMessageEdge(std::string filename,std::unordered_map<std::string,std::string> &person_0_0, 
                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_hasCreator_person_r, 
                    std::unordered_map<std::string,std::unordered_set<std::string>> &post_hasCreator_person_r,
                    std::unordered_map<std::string,std::string> &comment_hasCreator_person, 
                    std::unordered_map<std::string,std::string> &post_hasCreator_person,
                    std::unordered_map<std::string,std::unordered_set<std::string>> &post_hasTag_tag,
                    std::unordered_map<std::string,std::unordered_set<std::string>> &forum_containerOf_post_r,
                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_comment_r,
                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_post_r,
                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_comment,
                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_post) {
    std::string directory = filename;
    std::string line;

    std::cout<<"inputdir:"<<directory<<std::endl;

    // comment hasCreator person
    std::ifstream file(directory+"/comment_hasCreator_person_0_0.csv");
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }

        if (items.size() != 2) {
            std::cerr << "Invalid line: " << line << std::endl;
            continue;
        }
        comment_hasCreator_person_r[items[1]].insert(items[0]);
        comment_hasCreator_person[items[0]]=items[1];
    }
    file.close();

    // post hasCreator person
    file.open(directory+"/post_hasCreator_person_0_0.csv");
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }

        if (items.size() != 2) {
            std::cerr << "Invalid line: " << line << std::endl;
            continue;
        }
        post_hasCreator_person_r[items[1]].insert(items[0]);
        post_hasCreator_person[items[0]]=items[1];
    }
    file.close();

    // post hasTag tag
    file.open(directory+"/post_hasTag_tag_0_0.csv");
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }

        if (items.size() != 2) {
            std::cerr << "Invalid line: " << line << std::endl;
            continue;
        }
        post_hasTag_tag[items[0]].insert(items[1]);
    }
    file.close();

    file.open(directory+"/comment_replyOf_comment_0_0.csv");
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }

        if (items.size() != 2) {
            std::cerr << "Invalid line: " << line << std::endl;
            continue;
        }
        comment_replyOf_comment[items[0]].insert(items[1]);
        comment_replyOf_comment_r[items[1]].insert(items[0]);
        
    }

    file.close();

    file.open(directory+"/comment_replyOf_post_0_0.csv");
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }

        if (items.size() != 2) {
            std::cerr << "Invalid line: " << line << std::endl;
            continue;
        }
        comment_replyOf_post[items[0]].insert(items[1]);
        comment_replyOf_post_r[items[1]].insert(items[0]);
        
    }

    file.close();
    
	// forum_containerOf_post
    file.open(directory+"/forum_containerOf_post_0_0.csv");
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }

        if (items.size() != 2) {
            std::cerr << "Invalid line: " << line << std::endl;
            continue;
        }
        forum_containerOf_post_r[items[1]].insert(items[0]);
    }
    file.close();

	// person
    file.open(directory+"/person_0_0.csv");
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }
        person_0_0[items[0]] = line;
    }
    file.close();

    std::cout<<"load finished"<<std::endl;

}

void loadMessageCreationDate(std::string filename, std::unordered_map<std::string, std::string>& comment_creationDate, std::unordered_map<std::string, std::string>& post_creationDate) {
    std::string directory = filename;
    std::string line;
    std::ifstream file(directory+"/comment_0_0.csv");
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;
        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }
        comment_creationDate[items[0]] = items[1];
    }
    file.close();

    file.open(directory+"/post_0_0.csv");
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;
        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }
        post_creationDate[items[0]] = items[2];
    }
    file.close();
}

void loadMessageProperty(std::string filename, 
                 std::unordered_map<std::string, std::string>& comment_0_0, 
                 std::unordered_map<std::string, std::string>& post_0_0,
                 std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& person_likes_comment_r,
                 std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& person_likes_post_r,
                 std::unordered_map<std::string,std::unordered_set<std::string>>& post_isLocatedIn_place,
                 std::unordered_map<std::string,std::unordered_set<std::string>>& comment_hasTag_tag,
                 std::unordered_map<std::string,std::unordered_set<std::string>>& comment_isLocatedIn_place

                 ) {
    std::string directory = filename;
    std::string line;
	// comment
    std::ifstream file(directory+"/comment_0_0.csv");
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }
        comment_0_0[items[0]] = line;
    }
    file.close();

	// post
    file.open(directory+"/post_0_0.csv");
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }
        post_0_0[items[0]] = line;
    }
    file.close();

    // person_likes_comment
    file.open(directory+"/person_likes_comment_0_0.csv");
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }

        if (items.size() != 3) {
            std::cerr << "Invalid line: " << line << std::endl;
            continue;
        }
        person_likes_comment_r[items[1]][items[0]] = items[2];
    }
    file.close();

    // person_likes_post
    file.open(directory+"/person_likes_post_0_0.csv");
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }

        if (items.size() != 3) {
            std::cerr << "Invalid line: " << line << std::endl;
            continue;
        }
        person_likes_post_r[items[1]][items[0]] = items[2];
    }
    file.close();

	// post_isLocatedIn_place
    file.open(directory+"/post_isLocatedIn_place_0_0.csv");
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }

        if (items.size() != 2) {
            std::cerr << "Invalid line: " << line << std::endl;
            continue;
        }
        
        post_isLocatedIn_place[items[0]].insert(items[1]);
    }
    file.close();
	
	
	// comment_hasTag_tag
    file.open(directory+"/comment_hasTag_tag_0_0.csv");
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }

        if (items.size() != 2) {
            std::cerr << "Invalid line: " << line << std::endl;
            continue;
        }
        
        comment_hasTag_tag[items[0]].insert(items[1]);
    }
    file.close();
    
	// comment_isLocatedIn_place
    file.open(directory+"/comment_isLocatedIn_place_0_0.csv");
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }

        if (items.size() != 2) {
            std::cerr << "Invalid line: " << line << std::endl;
            continue;
        }
        comment_isLocatedIn_place[items[0]].insert(items[1]);
    }
    file.close();
  
}

// 初始化louvain算法输入带权图
void createLGraph(std::unordered_map<std::string, std::string> &person_0_0,
                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_hasCreator_person_r, 
                    std::unordered_map<std::string,std::unordered_set<std::string>> &post_hasCreator_person_r,
                    std::unordered_map<std::string,std::string> &comment_hasCreator_person, 
                    std::unordered_map<std::string,std::string> &post_hasCreator_person,
                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_comment,
                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_post,
                    Graph &g) {
    auto start = std::chrono::high_resolution_clock::now();
    // 处理每两个person点对
    for(auto i:person_0_0) {
        // 初始化person i的comment集
        std::unordered_set<std::string> iSet = comment_hasCreator_person_r[i.first]; 
        // dfs得到所有replyOf*邻居
        std::set<std::string> i_reply_comment,i_reply_post;
        for(auto ii: iSet) {
            std::unordered_map<std::string,std::unordered_map<std::string,int>> l_visited,l_visited1; 
            dfs_c(ii,comment_replyOf_comment,i_reply_comment,l_visited);
            dfs_c(ii,comment_replyOf_post,i_reply_post,l_visited1);
        }
        for(auto j: i_reply_comment) {
            // comment
            std::string to = comment_hasCreator_person[j];
            if(i.first!=to){
                g.addEdge(i.first,to,0.5);
            }
        }
        for(auto j: i_reply_post) {
            // post
            std::string to = post_hasCreator_person[j];
            if(i.first!=to){
                g.addEdge(i.first,to,0.5);
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    
    std::cout<<"create person reply graph finished,uses "<<duration.count()<<"ms"<<std::endl;
    
}

void divideCommunity(int part,long long messageNum,Graph& g,
                    std::unordered_map<std::string,std::string> &person_0_0, 
                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_hasCreator_person_r, 
                    std::unordered_map<std::string,std::unordered_set<std::string>> &post_hasCreator_person_r,
                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_comment,
                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_post,
                    unordered_map<string, std::unordered_set<std::string>>& communitiesSet,
                    std::vector<std::unordered_set<std::string>> &commentSetArray,
                    std::vector<std::unordered_set<std::string>> &postSetArray
                  ) {
    std::cout<<"divide community start"<<std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    std::unordered_set<std::string> communityMessageSet; //存储当前社区中person发布的message，用于判断阈值
    std::unordered_set<std::string> pPerson;  //存储已分配到社区中person
    int communityNo=0;
    int thresh = 0;
    while(g.getAdjList().size()>0) {
        // 取一个点加入社区
        auto it = g.getAdjList().begin();
        while(it!=g.getAdjList().end()&&pPerson.find(it->first)!=pPerson.end()) {
            it++;
        }
        if(it==g.getAdjList().end()) {
            break;
        }
        string center = it->first;
        pPerson.insert(center);
        communitiesSet[to_string(communityNo)].insert(center);
        communityMessageSet.insert(comment_hasCreator_person_r[center].begin(),comment_hasCreator_person_r[center].end());
        communityMessageSet.insert(post_hasCreator_person_r[center].begin(),post_hasCreator_person_r[center].end());
        thresh = communityMessageSet.size();
        if(thresh>=(messageNum/part)) {
            std::cout<<"divide community "<<communityNo<<" finished,community size:"<<communitiesSet[to_string(communityNo)].size()<<std::endl;
            communityNo++;
            thresh = 0;
            for(auto i:communitiesSet[to_string(communityNo)]) {
                g.removeNodes({i});
            }
            communityMessageSet.clear();
            continue;
        }
        PointHeap pointHeap;
        for(auto neighbor: g.getNeighbor(center)) {
            pointHeap.push({g.getNodeWeight(neighbor),neighbor});
        }
        // 取出社区中心的最大邻居，直到社区超过阈值或者没有邻居
        while(!pointHeap.isEmpty()) {
            // 取出最大邻居
            string tmp = pointHeap.pop().name;
            // 邻居已经分配到社区中就跳过
            if(pPerson.find(tmp)!=pPerson.end()) {
                continue;
            }
            pPerson.insert(tmp);
            communitiesSet[to_string(communityNo)].insert(tmp);
            communityMessageSet.insert(comment_hasCreator_person_r[tmp].begin(),comment_hasCreator_person_r[tmp].end());
            communityMessageSet.insert(post_hasCreator_person_r[tmp].begin(),post_hasCreator_person_r[tmp].end());
            thresh = communityMessageSet.size();
            if(thresh>=(messageNum/part)) {
                std::cout<<"divide community "<<communityNo<<" finished,community size:"<<communitiesSet[to_string(communityNo)].size()<<std::endl;
                communityNo++;
                thresh = 0;
                for(auto i:communitiesSet[to_string(communityNo)]) {
                    g.removeNodes({i});
                }
                communityMessageSet.clear();
                break;
            }
        }
    }
    // 分配不在社区中person，就是没有发布message的person或者message没有关联的person
    for(auto i:person_0_0) {
        if(pPerson.size()>=person_0_0.size()) {
            break;
        }
        if(pPerson.find(i.first)==pPerson.end()) {
            string cur = i.first;
            pPerson.insert(cur);
            communitiesSet[to_string(communityNo)].insert(cur);

            communityMessageSet.insert(comment_hasCreator_person_r[cur].begin(),comment_hasCreator_person_r[cur].end());
            communityMessageSet.insert(post_hasCreator_person_r[cur].begin(),post_hasCreator_person_r[cur].end());
            thresh = communityMessageSet.size();
            if(thresh>=(messageNum/part)) {
                std::cout<<"divide community "<<communityNo<<" finished,community size:"<<communitiesSet[to_string(communityNo)].size()<<std::endl;
                communityNo++;
                communityMessageSet.clear();
                thresh = 0;
            }
        }
    }
    if(communityMessageSet.size()>0)
        std::cout<<"divide community "<<communityNo<<" finished,community size:"<<communitiesSet[to_string(communityNo)].size()<<std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout<<"divide community finished,uses "<<duration.count()<<"s"<<endl;
    // 创建输出文件夹
    std::string outPath="./output";
    if (!createDirectory(outPath)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    outPath="./output/messages";
    if (!createDirectory(outPath)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    std::string routePath="./route";
    if (!createDirectory(routePath)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    std::ofstream message2partf("./route/message2part.csv");
    std::ofstream person2partf("./route/person2messagePart.csv");
    int commuCnt = 0;
    for(auto community:communitiesSet) {
        for(auto person:community.second) {
            int tmpCnt=0;
            person2partf<<person<<'|'<<community.first<<std::endl;
            for(auto i: comment_hasCreator_person_r[person]) {
                
                commentSetArray[std::stol(community.first)].insert(i);
                message2partf<<i<<'|'<<community.first<<std::endl;
            }
            for(auto i: post_hasCreator_person_r[person]) {
                
                postSetArray[std::stol(community.first)].insert(i);
                message2partf<<i<<'|'<<community.first<<std::endl;
            }
        }
        
    }
    message2partf.close();
    person2partf.close();
}

void divideMessageChain(unordered_map<string, std::unordered_set<std::string>>& communitiesSet,
                     std::vector<std::unordered_set<std::string>> &commentSetArray,
                     std::vector<std::unordered_set<std::string>> &postSetArray,
                     std::unordered_map<std::string,std::string> &comment_hasCreator_person, 
                     std::unordered_map<std::string,std::string> &post_hasCreator_person,
                     std::unordered_map<std::string,std::unordered_set<std::string>> &post_containerOf_forum,
                     std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_comment_r,
                     std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_post_r,
                     std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_comment,
                     std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_post,
                     std::unordered_map<std::string,std::unordered_set<std::string>> &post_hasTag_tag,
                     std::unordered_map<std::string,std::string> &comment_creationDate,
                     std::unordered_map<std::string,std::string> &post_creationDate) {
    std::string outPath="./output";
    for(int i=0;i<communitiesSet.size();i++){
            outPath="./output/messages/";
            outPath+=std::to_string(i);
            // 如果文件夹不存在，则创建它
            if (!createDirectory(outPath)) {
                // std::cerr << "Failed or no need to create directory."<<outPath << std::endl;
            }
            std::ofstream outfile_comment_hasCreator_person(outPath+"/comment_hasCreator_person_0_0.csv");
            std::ofstream outfile_post_hasCreator_person(outPath+"/post_hasCreator_person_0_0.csv");
            std::ofstream outfile_comment_0_0_added(outPath+"/comment_0_0_added.csv");
            std::ofstream outfile_post_0_0_added(outPath+"/post_0_0_added.csv");

            std::ofstream outfile_comment_replyOf_comment(outPath+"/comment_replyOf_comment_0_0.csv");
            std::ofstream outfile_comment_replyOf_post(outPath+"/comment_replyOf_post_0_0.csv");
            std::ofstream outfile_comment_replyOf_post_0_0_end(outPath+"/comment_replyOf_post_0_0_end.csv");

            std::ofstream outfile_post_hasTag_tag(outPath+"/post_hasTag_tag_0_0.csv");
            std::ofstream outfile_forum_containerOf_post(outPath+"/forum_containerOf_post_0_0.csv");
            std::ofstream outfile_forum(outPath+"/forum_0_0.csv");
            std::ofstream outfile_person(outPath+"/person_0_0.csv");
            // 写入person
            for(auto i:communitiesSet[std::to_string(i)]) {
                outfile_person<<i<<std::endl;
            }
            outfile_person.close();

            
            // 去重comment2comment
            std::unordered_map<std::string,std::unordered_map<std::string,int>> visited; 
            // 去重comment2post
            std::unordered_map<std::string,std::unordered_map<std::string,int>> visited1; 

            std::unordered_set<std::string> tmpComment; //当前分片的comment
            std::unordered_set<std::string> tmpPost;    //当前分片的post
            std::unordered_set<std::string> tmpForum;   //当前分片的forum
            // 存储commentPath (fromComment neighborComment endComment)
            std::unordered_map<std::string,std::vector<std::string>> commentPath;

            tmpComment.insert(commentSetArray[i].begin(),commentSetArray[i].end());
            tmpPost.insert(postSetArray[i].begin(),postSetArray[i].end());

            // 开始划分message
            for(auto k:commentSetArray[i]) {
                // comment的前向comment和边
                tmpComment.insert(comment_replyOf_comment_r[k].begin(),comment_replyOf_comment_r[k].end());
                for(auto kk:comment_replyOf_comment_r[k]){
                    if(visited[kk].find(k)==visited[kk].end()) {
                            visited[kk][k] = 1;
                            outfile_comment_replyOf_comment<<kk<<'|'<<k<<std::endl;
                    }
                }
                // comment的后向comment和边
                tmpComment.insert(comment_replyOf_comment[k].begin(),comment_replyOf_comment[k].end());
                for(auto kk:comment_replyOf_comment[k]){
                    if(visited[k].find(kk)==visited[k].end()) {
                            visited[k][kk] = 1;
                            outfile_comment_replyOf_comment<<k<<'|'<<kk<<std::endl;
                    }
                }
                // comment的后向post和边
                tmpPost.insert(comment_replyOf_post[k].begin(),comment_replyOf_post[k].end());
                for(auto kk:comment_replyOf_post[k]){
                    if(visited1[k].find(kk)==visited1[k].end()) {
                            visited1[k][kk] = 1;
                            outfile_comment_replyOf_post<<k<<'|'<<kk<<std::endl;
                    }
                }
            }
            for(auto k:postSetArray[i]) {
                // post的前向
                tmpComment.insert(comment_replyOf_post_r[k].begin(),comment_replyOf_post_r[k].end());
                for(auto kk:comment_replyOf_post_r[k]){
                    if(visited1[kk].find(k)==visited1[kk].end()) {
                            visited1[kk][k] = 1;
                            outfile_comment_replyOf_post<<kk<<'|'<<k<<std::endl;
                    }
                }
            }
            visited.clear();
            visited1.clear();
            
            for (auto k:commentSetArray[i]) {
                std::vector<std::string> tmp;
                dfs(k, comment_replyOf_comment, tmp, commentPath);
            }
            for(auto k=commentPath.begin();k!=commentPath.end();k++){
                std::string from = k->first;
                std::vector<std::string> tmp = k->second;
                if(tmp.size()==1) {
                    // comment没有comment邻居
                    tmpComment.insert(tmp[0]);  // useless
                    for(std::string j : comment_replyOf_post[tmp[0]]) {
                        tmpPost.insert(j);
                        if(visited1[tmp[0]].find(j)==visited1[tmp[0]].end()) {
                            visited1[tmp[0]][j] = 1;
                            outfile_comment_replyOf_post_0_0_end<<tmp[0]<<'|'<<j<<std::endl;
                        }
                    }
                } else if(tmp.size()==2) {
                    // comment只有一步comment
                    tmpComment.insert(tmp[1]);
                    if(visited[tmp[0]].find(tmp[1])==visited[tmp[0]].end()) {
                        visited[tmp[0]][tmp[1]] = 1;
                        outfile_comment_replyOf_comment<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                    }
                    for(std::string j : comment_replyOf_post[tmp[1]]) {
                        // comment存在comment链
                        tmpPost.insert(j);
                        if(visited1[tmp[0]].find(j)==visited1[tmp[0]].end()) {
                            visited1[tmp[0]][j] = 1;
                            outfile_comment_replyOf_post_0_0_end<<tmp[0]<<'|'<<j<<std::endl;
                        }
                    }
                } else {
                    // 顶点有多步邻居
                    tmpComment.insert(tmp[1]);
                    if(visited[tmp[0]].find(tmp[1])==visited[tmp[0]].end()) {
                        visited[tmp[0]][tmp[1]] = 1;
                        outfile_comment_replyOf_comment<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                    }
                    // 这里把末端邻居合并为post
                    for(std::string j : comment_replyOf_post[tmp[2]]) {
                        tmpPost.insert(j);
                        if(visited1[tmp[0]].find(j)==visited1[tmp[0]].end()) {
                            visited1[tmp[0]][j] = 1;
                            outfile_comment_replyOf_post_0_0_end<<tmp[0]<<'|'<<j<<std::endl;
                        }
                        // 压缩了路径，不过会导致comment2post变多
                    }
                }
            }

            // 写入comment
            for(auto k:tmpComment){
                if(commentSetArray[i].find(k)!=commentSetArray[i].end()) {
                    outfile_comment_hasCreator_person<<k<<'|'<<comment_hasCreator_person[k]<<std::endl;
                } else {
                    outfile_comment_0_0_added<<k<<"|"<<comment_creationDate[k]<<std::endl;
                }
            }
            // 写入post
            for(auto k:tmpPost){
                if(postSetArray[i].find(k)!=postSetArray[i].end()) {
                    outfile_post_hasCreator_person<<k<<'|'<<post_hasCreator_person[k]<<std::endl;
                } else {
                    outfile_post_0_0_added<<k<<"|"<<post_creationDate[k]<<std::endl;
                }
                for(auto j:post_hasTag_tag[k]){
                    outfile_post_hasTag_tag<<k<<'|'<<j<<std::endl;
                }
                for(auto j:post_containerOf_forum[k]){
                    outfile_forum_containerOf_post<<j<<'|'<<k<<std::endl;
                    tmpForum.insert(j);
                }
            }
            
            // 写入forum
            for(auto k:tmpForum){
                outfile_forum<<k<<std::endl;
            }

            // 释放空间
            visited.clear();
            visited.rehash(0);
            visited1.clear();
            visited1.rehash(0);
            tmpComment.clear();
            tmpComment = std::unordered_set<string>();
            tmpPost.clear();
            tmpPost = std::unordered_set<string>();
            tmpForum.clear();
            tmpForum = std::unordered_set<string>();
            for (auto& pair : commentPath) {
                pair.second.clear();
                pair.second.shrink_to_fit();
            }
            commentPath.clear();
            commentPath.rehash(0);
            
            outfile_comment_hasCreator_person.close();
            outfile_post_hasCreator_person.close();
            outfile_comment_0_0_added.close();
            outfile_post_0_0_added.close();
            outfile_comment_replyOf_comment.close();
            outfile_comment_replyOf_post.close();
            outfile_post_hasTag_tag.close();
            outfile_forum_containerOf_post.close();
            outfile_forum.close();
        }
}

void divideMessageProperty(std::vector<std::unordered_set<std::string>> &commentSetArray,
                   std::vector<std::unordered_set<std::string>> &postSetArray,
                   std::unordered_map<std::string,std::string> &comment_0_0,
                   std::unordered_map<std::string,std::string> &post_0_0,
                   std::unordered_map<std::string,std::unordered_set<std::string>> &comment_hasTag_tag,
                   std::unordered_map<std::string,std::unordered_set<std::string>> &comment_isLocatedIn_place,
                   std::unordered_map<std::string,std::unordered_set<std::string>> &post_isLocatedIn_place,
                   std::unordered_map<std::string,std::unordered_map<std::string,std::string>> &person_likes_comment_r,
                   std::unordered_map<std::string,std::unordered_map<std::string,std::string>> &person_likes_post_r) {
    string outPath = "./output/";
    for(int i=0;i<commentSetArray.size();i++){
        outPath="./output/messages/";
        outPath+=std::to_string(i);
        // 如果文件夹不存在，则创建它
        if (!createDirectory(outPath)) {
            // std::cerr << "Failed or no need to create directory."<<outPath << std::endl;
        }
        std::ofstream outfile_comment(outPath+"/comment_0_0.csv");
        std::ofstream outfile_post(outPath+"/post_0_0.csv");
        
        std::ofstream outfile_comment_isLocatedIn_place(outPath+"/comment_isLocatedIn_place_0_0.csv");
        std::ofstream outfile_post_isLocatedIn_place(outPath+"/post_isLocatedIn_place_0_0.csv");
        std::ofstream outfile_comment_hasTag_tag(outPath+"/comment_hasTag_tag_0_0.csv");

        std::ofstream outfile_person_likes_comment(outPath+"/person_likes_comment_0_0.csv");
        std::ofstream outfile_person_likes_post(outPath+"/person_likes_post_0_0.csv");
        

        // 开始划分message
        for(auto k:commentSetArray[i]) {
            outfile_comment<<comment_0_0[k]<<std::endl;
            for(auto j:comment_isLocatedIn_place[k]){
                outfile_comment_isLocatedIn_place<<k<<'|'<<j<<std::endl;
            }
            for(auto j:comment_hasTag_tag[k]){
                outfile_comment_hasTag_tag<<k<<'|'<<j<<std::endl;
            }
            for(auto j:person_likes_comment_r[k]){
                outfile_person_likes_comment<<j.first<<'|'<<k<<'|'<<j.second<<std::endl;
            }
        }
        for(auto k:postSetArray[i]) {
            outfile_post<<post_0_0[k]<<std::endl;
            for(auto j:post_isLocatedIn_place[k]){
                outfile_post_isLocatedIn_place<<k<<'|'<<j<<std::endl;
            }
            for(auto j:person_likes_post_r[k]){
                outfile_person_likes_post<<j.first<<'|'<<k<<'|'<<j.second<<std::endl;
            }
        }
        
        outfile_comment.close();
        outfile_post.close();
        outfile_comment_isLocatedIn_place.close();
        outfile_post_isLocatedIn_place.close();
        outfile_comment_hasTag_tag.close();

        outfile_person_likes_comment.close();
        outfile_person_likes_post.close();
        
    }
}

void forumAdd(std::unordered_map<std::string, std::string> &post_0_0) {
    std::string outPath="./output/forum";
    if (!createDirectory(outPath)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    std::ofstream outfile_post(outPath+"/post_0_0.csv");
    for(auto i:post_0_0) {
        outfile_post<<i.first<<std::endl;
    }
    outfile_post.close();
}

// 释放容器空间的函数
void releaseMemory(std::unordered_map<std::string, std::string>& map) {
    map.clear();
    std::unordered_map<std::string, std::string>().swap(map);
}

// 释放嵌套的 unordered_map 的空间
void releaseNestedMap(std::unordered_map<std::string, std::unordered_map<std::string, std::string>>& map) {
    for (auto& outer_pair : map) {
        outer_pair.second.clear();
        std::unordered_map<std::string, std::string>().swap(outer_pair.second);
    }
    map.clear();
    std::unordered_map<std::string, std::unordered_map<std::string, std::string>>().swap(map);
}

// 释放嵌套的 unordered_map 的空间
void releaseNestedMap(std::unordered_map<std::string, std::unordered_set<std::string>>& map) {
    for (auto& outer_pair : map) {
        outer_pair.second.clear();
        std::unordered_set<std::string>().swap(outer_pair.second);
    }
    map.clear();
    std::unordered_map<std::string, std::unordered_set<std::string>>().swap(map);
}


int main(int argc, char *argv[]) {
    
    auto start = std::chrono::high_resolution_clock::now();
    // 帮助信息
    std::string helpMessage = R"(
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
  ./program_name 0.1                       # Uses SCALE_FACTOR = 0.1, default DIRECTORY and PART
  ./program_name 0.1 /path/to/dir          # Uses SCALE_FACTOR = 0.1, DIRECTORY = /path/to/dir, and default PART
  ./program_name 0.1 /path/to/dir 5        # Uses SCALE_FACTOR = 0.1, DIRECTORY = /path/to/dir, and PART = 5
)";

    // 默认值
    std::string initialDir = "/home/shared/data/social_network-csv_composite-longdateformatter-sf0.1/dynamic/";
    std::string sfvalue = "0.1";
    int part = 3;

    // 检查是否请求帮助
    if (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
        std::cout << helpMessage << std::endl;
        return 0;
    }

    // 根据输入参数调整值
    if (argc >= 3) {
        sfvalue = argv[1];
        initialDir = argv[2];
    }

    if (argc == 4) {
        part = std::atoi(argv[3]);  // 将字符串转换为整数
    }
    
    if(argc > 4) {
        std::cout << helpMessage << std::endl;
        return 0;
    }

    // 输出选择的值（仅作为示例，实际代码中可以执行其他操作）
    std::cout << "Scale Factor: " << sfvalue << std::endl;
    std::cout << "Initial Directory: " << initialDir << std::endl;
    std::cout << "Partition Part Num: " << part << std::endl;

    // 存储数据
    std::unordered_map<std::string,std::unordered_set<std::string>> person_hasInterest_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> person_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_studyAt_organisation;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_workAt_organisation;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_likes_comment_r;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_likes_post_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasCreator_person_r;
    std::unordered_map<std::string,std::string> comment_hasCreator_person;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasCreator_person_r;
    std::unordered_map<std::string,std::string> post_hasCreator_person;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasTag_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasTag_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_set<std::string>> forum_containerOf_post;
    std::unordered_map<std::string,std::unordered_set<std::string>> forum_containerOf_post_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_replyOf_comment_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_replyOf_post_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_replyOf_comment;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_replyOf_post;
    std::unordered_map<std::string, std::string> person_0_0;
    std::unordered_map<std::string, std::string> comment_0_0;
    std::unordered_map<std::string, std::string> post_0_0;
    std::unordered_map<std::string, std::string> comment_creationDate;
    std::unordered_map<std::string, std::string> post_creationDate;
    std::unordered_map<std::string, std::string> forum_0_0;

    Graph g;    //存储person-replyOf带权图
    unordered_map<string, std::unordered_set<std::string>> communitiesSet;
    std::vector<std::unordered_set<std::string>> commentSetArray(part);
    std::vector<std::unordered_set<std::string>> postSetArray(part);
    std::map<std::string,long long> sf2messageNum = {
        {"0.1", 286744},
        {"3", 9010236},
        {"30", 87095182},
        {"100", 331896459},
        {"300", 956119240},
        {"1000", 3126402941}
    };
    // 获取划分开始时间点
    auto partitionStart = std::chrono::high_resolution_clock::now();
    // step1 划分person社区和messsage分片拓扑
    loadMessageEdge(initialDir,person_0_0,comment_hasCreator_person_r,post_hasCreator_person_r,comment_hasCreator_person,post_hasCreator_person,post_hasTag_tag,forum_containerOf_post_r,comment_replyOf_comment_r,comment_replyOf_post_r,comment_replyOf_comment,comment_replyOf_post);
    loadMessageCreationDate(initialDir,comment_creationDate,post_creationDate);
    createLGraph(person_0_0,comment_hasCreator_person_r,post_hasCreator_person_r,comment_hasCreator_person,post_hasCreator_person,comment_replyOf_comment,comment_replyOf_post,g);
    divideCommunity(part, sf2messageNum[sfvalue], g, person_0_0, comment_hasCreator_person_r, post_hasCreator_person_r,comment_replyOf_comment,comment_replyOf_post,communitiesSet,commentSetArray,postSetArray);
    divideMessageChain(communitiesSet,commentSetArray,postSetArray,comment_hasCreator_person,post_hasCreator_person,forum_containerOf_post_r,comment_replyOf_comment_r,comment_replyOf_post_r,comment_replyOf_comment,comment_replyOf_post,post_hasTag_tag,comment_creationDate,post_creationDate);
    // 释放空间
    releaseMemory(person_0_0);
    releaseNestedMap(comment_hasCreator_person_r);
    releaseNestedMap(post_hasCreator_person_r);
    releaseMemory(comment_hasCreator_person);
    releaseMemory(post_hasCreator_person);
    releaseNestedMap(comment_replyOf_comment);
    releaseNestedMap(comment_replyOf_post);
    releaseNestedMap(comment_replyOf_comment_r);
    releaseNestedMap(comment_replyOf_post_r);
    releaseNestedMap(forum_containerOf_post_r);
    releaseNestedMap(communitiesSet);
    auto firstStepEnd = std::chrono::high_resolution_clock::now();  //第一阶段结束

    auto firstDuration = std::chrono::duration_cast<std::chrono::seconds>(firstStepEnd - partitionStart);
    std::cout << "partition step1 finished,uses " << firstDuration.count() << " seconds" << std::endl;

    // step2 中心Message属性和Likes
    loadMessageProperty(initialDir, comment_0_0, post_0_0,person_likes_comment_r,person_likes_post_r,post_isLocatedIn_place,comment_hasTag_tag,comment_isLocatedIn_place); 
    divideMessageProperty(commentSetArray,postSetArray,comment_0_0,post_0_0,comment_hasTag_tag,comment_isLocatedIn_place,post_isLocatedIn_place,person_likes_comment_r,person_likes_post_r);
    // forum分片
    forumAdd(post_0_0);

    // 获取结束时间点
    auto end = std::chrono::high_resolution_clock::now();
    auto secondDuration = std::chrono::duration_cast<std::chrono::seconds>(end - firstStepEnd);
    std::cout << "partition step2 finished,uses " << secondDuration.count() << " seconds" << std::endl;


    // 计算所消耗的时间
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << "Total time: " << duration.count() << " seconds" << std::endl;


    return 0;
}