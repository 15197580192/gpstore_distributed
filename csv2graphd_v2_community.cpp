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
// std::vector<std::thread> threads; // 用于存储线程的向量，其实没用

void printMemoryUsage() {
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        double memoryUsageMB = usage.ru_maxrss / 1024.0; // 将 KB 转换为 MB
        std::cout << "Memory Usage: " << memoryUsageMB << " MB" << std::endl;
    } else {
        std::cerr << "Failed to get memory usage" << std::endl;
    }
}

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


class Person {
public:
    std::string id;
    std::string line;

    void printInfo() const {
        std::cout << "ID: " << id << std::endl;
        std::cout << "line: " << line << std::endl;
        std::cout << "------------------------" << std::endl;
    }
};


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
            // if(kv.second.empty()) {
            //     adjList[kv.first].push_back({kv.first, kv.first, 0.0});
            // }
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

// Class to represent Disjoint Set
class DisjointSet {
    std::map<std::string, std::string> parent;

public:
    // Constructor to initialize parent of each element to itself
    DisjointSet(const std::set<std::string>& elements) {
        for (const std::string& e : elements) {
            parent[e] = e;
        }
    }

    // Find operation to find the root of the element
    std::string find(const std::string& x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);
        }
        return parent[x];
    }

    // Union operation to merge two sets
    void merge(const std::string& x, const std::string& y) {
        std::string rootX = find(x);
        std::string rootY = find(y);
        if (rootX != rootY) {
            parent[rootY] = rootX;
        }
    }
};

// 定义线程安全队列
class SafeQueue {
public:
    void push(const std::string& value) {
        std::lock_guard<std::mutex> lock(mtx);
        q.push(value);
        cv.notify_one();
    }

    bool pop(std::string& value) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&]() { return !q.empty() || done; });
        if (!q.empty()) {
            value = q.front();
            q.pop();
            return true;
        }
        return false;
    }

    void setDone() {
        std::lock_guard<std::mutex> lock(mtx);
        done = true;
        cv.notify_all();
    }

private:
    std::queue<std::string> q;
    std::mutex mtx;
    std::condition_variable cv;
    bool done = false;
};

// 将数据写入文件
void writer(SafeQueue& queue, const std::string& filename) {
    // 创建一个文件输出流
    std::ofstream file(filename, std::ios::out | std::ios::binary);
    // 如果文件无法打开，输出错误信息并返回
    if (!file) {
        std::cerr << "File could not be opened!" << std::endl;
        return;
    }

    std::string str;
    while (queue.pop(str)) {
        file.write(str.c_str(), str.size());
    }

    // 关闭文件
    file.close();
}

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

void getKnows(std::string givenPersonID, std::map<std::string, std::map<std::string,std::string>>& knows, std::set<std::string>& twoHopResults) {
    // Get two-hop knows relationships for a given person ID
    twoHopResults.insert(givenPersonID);

    // 获取文件夹路径
    // std::string directory = filename.substr(0, filename.find_last_of('/'));

    // 如果文件夹不存在，则创建它
    // if (!createDirectory(directory)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    // }

    // std::ofstream ofile(filename);
    // ofile << "Person.id|Person.id|creationDate" << std::endl;
    // std::cout<<knows.size();

    for (const auto& edge : knows[givenPersonID]) {
        std::string target = edge.first;
        std::string content = knows[givenPersonID][target];
        twoHopResults.insert(target);
        // ofile <<givenPersonID<<'|'<<target<<'|'<< content<<std::endl;

        for (const auto& edge2 : knows[target]) {
            std::string target2 = edge2.first;
            content = knows[target][target2];
            twoHopResults.insert(target2);
            // ofile <<target<<'|'<<target2<<'|'<< content<<std::endl;
        }
    }
    // ofile.close();
}

void getPerson(std::set<std::string> twoHopResults, std::string filename, std::map<std::string, Person> people) {
    // 获取文件夹路径
    std::string directory = filename.substr(0, filename.find_last_of('/'));

    // 如果文件夹不存在，则创建它
    if (!createDirectory(directory)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }

    std::ofstream ofile(filename);
    ofile << "Person.id|Person.firstName|Person.lastName|Person.birthday|Person.creationDate|Person.gender|Person.browserUsed|Person.locationIP" << std::endl;
    for(auto it=twoHopResults.begin();it!=twoHopResults.end();it++){
        ofile<<people[*it].line<<std::endl;
    }
}

void dfs(std::string current, std::unordered_map<std::string, std::unordered_set<std::string>> &reply1, std::vector<std::string> &tmp, std::unordered_map<std::string,std::vector<std::string>> &commentPath) {
//    需要 from neighbor end post
    if(tmp.size()<=2) {
        // 最多存三个，用于合并后续邻居(起点，一步邻居，末端comment)
        tmp.push_back(current);
    }
    if(reply1[current].size()==0) {
        if(tmp.size()==1) {
            // 顶点没有邻居
        } else if(tmp.size()==2) {
            // 顶点只有一步邻居
        } else {
            // if(tmp[2]!=current)
            //     std::cout<<tmp[0]<<' '<<tmp[1]<<' '<<tmp[2]<<' '<<current<<std::endl;
            // 顶点有多步邻居
            tmp.pop_back();
            tmp.push_back(current);
        }
        commentPath[tmp[0]] = tmp;
        return;
    }
    for (auto i : reply1[current]) {
        dfs(i, reply1, tmp, commentPath);
    }
}

void dfs1(std::string outfile, std::string current, std::unordered_map<std::string, std::unordered_set<std::string>> &reply1, std::unordered_set<std::string> &end) {
//    需要 from neighbor end post
    if(reply1[current].size()==0) {
        end.insert(current);
        return;
    }
    for (auto i : reply1[current]) {
        std::ofstream ofile(outfile, std::ios::app);
        ofile<<i<<'|'<<current<<std::endl;  
        ofile.close();
        dfs1(outfile, i, reply1, end);
    }
}

// 检查边是否已访问
bool isEdgeVisited(const std::string& from, const std::string& to, const std::map<std::string, std::map<std::string,int>>& visited) {
    auto fromIt = visited.find(from);
    if (fromIt != visited.end()) {
        auto toIt = fromIt->second.find(to);
        if (toIt != fromIt->second.end()) {
            return toIt->second == 1;
        }
    }
    return false;
}

void getDivide(std::set<std::string> &twoHopResults, std::string filename, std::string outdir, std::set<std::string>& comment, std::set<std::string>& post, std::map<std::string, std::map<std::string,std::string>>& knows) {
    // 获取文件夹路径
    // std::string directory = filename.substr(0, filename.find_last_of('/'));
    std::string directory = filename;
    std::string directory1 = "./0/";
    directory1 = outdir;

    // 如果文件夹不存在，则创建它
    if (!createDirectory(directory1)) {
        std::cerr << "Failed or no need to create directory." << std::endl;
    }
    std::cout<<"inputdir:"<<directory<<std::endl;
    std::cout<<"outputdir:"<<directory1<<std::endl;
    
    std::unordered_map<std::string, std::unordered_set<std::string>> reply1;
    std::unordered_map<std::string, std::unordered_set<std::string>> reply2;

    // std::unordered_map<std::string, std::unordered_map<std::string, std::string>> creator1;
    // std::unordered_map<std::string, std::unordered_map<std::string, std::string>> creator2;
    
    
    std::set<std::string> forum;
    
    // person_hasInterest_tag
    
    std::ifstream file(directory+"/person_hasInterest_tag_0_0.csv");
    std::ofstream file2(directory1+"/person_hasInterest_tag_0_0.csv");
    
    std::string line;

    // std::getline(file, line); // Skip header line
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
        
        if(twoHopResults.find(items[0])!=twoHopResults.end()){
            file2<<line<<std::endl;
        }
    }
    file.close();
    file2.close(); 
    
    // person_isLocatedIn_city
    
    file.open(directory+"/person_isLocatedIn_place_0_0.csv");
    file2.open(directory1+"/person_isLocatedIn_place_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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
        
        if(twoHopResults.find(items[0])!=twoHopResults.end()){
            file2<<line<<std::endl;
        }
    }
    file.close();
    file2.close(); 
    
    // person_studyAt_university
    
    file.open(directory+"/person_studyAt_organisation_0_0.csv");
    file2.open(directory1+"/person_studyAt_organisation_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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
        
        if(twoHopResults.find(items[0])!=twoHopResults.end()){
            file2<<line<<std::endl;
        }
    }
    file.close();
    file2.close(); 
    
    
    
    // person_workAt_company
    
    
    file.open(directory+"/person_workAt_organisation_0_0.csv");
    file2.open(directory1+"/person_workAt_organisation_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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
        
        if(twoHopResults.find(items[0])!=twoHopResults.end()){
            file2<<line<<std::endl;
        }
    }
    file.close();
    file2.close(); 
    
    
    // forum_hasMember_person
    
    file.open(directory+"/forum_hasMember_person_0_0.csv");
    file2.open(directory1+"/forum_hasMember_person_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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
        
        if(twoHopResults.find(items[1])!=twoHopResults.end()){
            file2<<line<<std::endl;
            forum.insert(items[0]);
        }
    }
    file.close();
    file2.close(); 
    
    // forum_hasModerator_person

    file.open(directory+"/forum_hasModerator_person_0_0.csv");
    file2.open(directory1+"/forum_hasModerator_person_0_0.csv");
    std::map<std::string,std::set<std::string>> forum2person;
    
    // std::getline(file, line); // Skip header line
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
        forum2person[items[0]].insert(items[1]);
        
        if(twoHopResults.find(items[1])!=twoHopResults.end()){
            file2<<line<<std::endl;
            forum.insert(items[0]);
        }
    }
    file.close();
    file2.close(); 

    // comment hasCreator person
    file.open(directory+"/comment_hasCreator_person_0_0.csv");
    // 这里不用输出，后面统一输出
    // file2.open(directory1+"/comment_hasCreator_person_0_0.csv");
    // file2<<"Comment.id|Person.id"<<std::endl;    //为什么没有输出到呢，open会覆盖

    // std::getline(file, line); // Skip header line
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
        
        if(twoHopResults.find(items[1])!=twoHopResults.end()){
            comment.insert(items[0]);
            // file2<<line<<std::endl;
        }
    }
    file.close();
    // file2.close();

    // post hasCreator person
    file.open(directory+"/post_hasCreator_person_0_0.csv");
    file2.open(directory1+"/post_hasCreator_person_0_0.csv");
    // file2<<"Post.id|Person.id"<<std::endl;

    // std::getline(file, line); // Skip header line
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
        
        if(twoHopResults.find(items[1])!=twoHopResults.end()){
            post.insert(items[0]);
            file2<<line<<std::endl;
        }
    }
    file.close();
    file2.close();

    file.open(directory+"/comment_replyOf_comment_0_0.csv");
    // std::getline(file, line); // Skip header line
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

        reply1[items[0]].insert(items[1]);
    }
    file.close();

    file.open(directory+"/comment_replyOf_post_0_0.csv");
    // std::getline(file, line); // Skip header line
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

        reply2[items[0]].insert(items[1]);
    }
    file.close();

    // // version1:commet链至post
    // // bfs or dfs?
    // // 无需考虑新引入的commet? todo
    // // file2.open(directory1+"comment_replyOf_comment_0_0.csv");
    // std::ofstream file3(directory1+"comment_replyOf_post_0_0.csv");
    // // file3<<"Comment.id|Post.id"<<std::endl;
    // std::unordered_set<std::string> end;
    // // 所有的comment commment-neigbor post
    // for (auto i=comment.begin();i!=comment.end();i++) {
    //     std::vector<std::string> tmp;
    //     dfs1(directory1+"comment_replyOf_comment_0_0.csv", *i, reply1, end);
	// }
    // for(auto i=end.begin();i!=end.end();i++) {
    // 	for(std::string j : reply2[*i]) {
    //         file3<<*i<<'|'<<j<<std::endl;
    //     }
    // }
    // file3.close();

	// version2:commet链至post
    // bfs or dfs?
    // 无需考虑新引入的commet? todo
    // std::ofstream file3(directory1+"comment_replyOf_post_0_0.csv");
    std::ofstream file3_1(directory1+"/comment_replyOf_comment_0_0.csv");
    std::ofstream file3_2(directory1+"/comment_replyOf_post_0_0.csv");
    // 去重comment2comment
    std::unordered_map<std::string,std::unordered_map<std::string,int>> visited; 
    // 去重comment2post
    std::unordered_map<std::string,std::unordered_map<std::string,int>> visited1; 
    // file3<<"Comment.id|Post.id"<<std::endl;
    std::unordered_set<std::string> end;
    // 所有的comment commment-neigbor post
    std::unordered_map<std::string,std::vector<std::string>> commentPath;
    for (auto i=comment.begin();i!=comment.end();i++) {
        std::vector<std::string> tmp;
        dfs(*i, reply1, tmp, commentPath);
	}
    // std::set<std::string> comment1=comment;
    // std::set<std::string> post1=post;

    for(auto i=commentPath.begin();i!=commentPath.end();i++) {
        std::string from = i->first;
    	std::vector<std::string> tmp = i->second;
        if(tmp.size()==1) {
            // comment没有comment邻居
            // ！！！这里也要存post
            comment.insert(tmp[0]); // 里面本来就存了
            // file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
            for(std::string j : reply2[tmp[0]]) {
                post.insert(j);
                if(visited1[tmp[0]].find(j)==visited1[tmp[0]].end()) {
                    visited1[tmp[0]][j] = 1;
                    file3_2<<tmp[0]<<'|'<<j<<std::endl;
                }
                // file3_2<<tmp[0]<<'|'<<j<<std::endl;
            }
        } else if(tmp.size()==2) {
            // comment只有一步comment
            comment.insert(tmp[1]);
            if(visited[tmp[0]].find(tmp[1])==visited[tmp[0]].end()) {
                visited[tmp[0]][tmp[1]] = 1;
                file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
            }
            // file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
            for(std::string j : reply2[tmp[1]]) {
                // comment存在comment链
                post.insert(j);
                if(visited1[tmp[1]].find(j)==visited1[tmp[1]].end()) {
                    visited1[tmp[1]][j] = 1;
                    file3_2<<tmp[1]<<'|'<<j<<std::endl;
                }
                // file3_2<<tmp[1]<<'|'<<j<<std::endl;
            }
        } else {
            // 顶点有多步邻居
            comment.insert(tmp[1]);
            if(visited[tmp[0]].find(tmp[1])==visited[tmp[0]].end()) {
                visited[tmp[0]][tmp[1]] = 1;
                file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
            }
            // file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
            // 这里把末端邻居合并为post
            for(std::string j : reply2[tmp[2]]) {
                post.insert(j);
                if(visited1[tmp[1]].find(j)==visited1[tmp[1]].end()) {
                    visited1[tmp[1]][j] = 1;
                    file3_2<<tmp[1]<<'|'<<j<<std::endl;
                }
                // 压缩了路径，不过好像导致comment2post变多了？
            }
        }
    }
	
	// comment_hasCreator_person
    file.open(directory+"/comment_hasCreator_person_0_0.csv");
    file2.open(directory1+"/comment_hasCreator_person_0_0.csv");// 这里加了所有comment到person的边，所以之前不用输出边，这里覆盖是对的
    // file2<<"Comment.id|Person.id"<<std::endl;
    std::set <std::string> newPerson;
    
    // std::getline(file, line); // Skip header line
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
        
        if(comment.find(items[0])!=comment.end()){
        	if(twoHopResults.find(items[1])==twoHopResults.end()) {
        		
           		twoHopResults.insert(items[1]);
           		newPerson.insert(items[1]);
			}
            file2<<line<<std::endl;
        }
    }
    file.close();
    file2.close();
    
    
	// post_hasCreator_person
    file.open(directory+"/post_hasCreator_person_0_0.csv");
    file2.open(directory1+"/post_hasCreator_person_0_0.csv", std::ios::app);
    // std::set <std::string> newPerson;
    
    // std::getline(file, line); // Skip header line
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
        
        if(post.find(items[0])!=post.end()){
        	if(twoHopResults.find(items[1])==twoHopResults.end()) {
        		
           		twoHopResults.insert(items[1]);
           		newPerson.insert(items[1]);
			}
            file2<<line<<std::endl;
        }
    }
    file.close();
    file2.close();
    
    // todo:newPerson与其他person的knows,finished!
    // 是否需要传入knows? 传入吧 直接打开新建也是一样
    // 遍历twoHopResults中的每个person
    file2.open(directory1+"/person_knows_person_0_0.csv");
    std::map<std::string, std::map<std::string,int>> kvisited;
    for (const auto& person : twoHopResults) {
        if (knows.find(person) != knows.end()) {
            // 找到person的所有邻居
            for (const auto& neighbor : knows[person]) {
                // 将每个person与其邻居之间的knows边存储起来
                if(twoHopResults.find(neighbor.first)!=twoHopResults.end()){
                    if(!isEdgeVisited(neighbor.first, person, kvisited)) {
                        kvisited[person][neighbor.first] = 1;
                        file2<<person<<'|'<<neighbor.first<<'|'<<neighbor.second<<std::endl;
                    }
                }
            }
        }
    }
    file2.close();
    
    
	// post_hasTag_tag
    file.open(directory+"/post_hasTag_tag_0_0.csv");
    file2.open(directory1+"/post_hasTag_tag_0_0.csv");
    // std::set <std::string> newPerson;
    
    // std::getline(file, line); // Skip header line
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
        
        if(post.find(items[0])!=post.end()){
            file2<<line<<std::endl;
        }
    }
    file.close();
    file2.close(); 

	// post_isLocatedIn_place
    file.open(directory+"/post_isLocatedIn_place_0_0.csv");
    file2.open(directory1+"/post_isLocatedIn_place_0_0.csv");
    // std::set <std::string> newPerson;
    
    // std::getline(file, line); // Skip header line
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
        
        if(post.find(items[0])!=post.end()){
            file2<<line<<std::endl;
        }
    }
    file.close();
    file2.close(); 
	
	
	// comment_hasTag_tag
    file.open(directory+"/comment_hasTag_tag_0_0.csv");
    file2.open(directory1+"/comment_hasTag_tag_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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
        
        if(comment.find(items[0])!=comment.end()){
            file2<<line<<std::endl;
        }
    }
    file.close();
    file2.close(); 
    
	// comment_isLocatedIn_place
    file.open(directory+"/comment_isLocatedIn_place_0_0.csv");
    file2.open(directory1+"/comment_isLocatedIn_place_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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
        
        if(comment.find(items[0])!=comment.end()){
            file2<<line<<std::endl;
        }
    }
    file.close();
    file2.close(); 
    
    
	// forum_containerOf_post
    file.open(directory+"/forum_containerOf_post_0_0.csv");
    file2.open(directory1+"/forum_containerOf_post_0_0.csv");
    std::set <std::string> newForum;
    
    // std::getline(file, line); // Skip header line
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
        
        if(post.find(items[1])!=post.end()){
            file2<<line<<std::endl;
            forum.insert(items[0]); // 提前已经完成了person_forum,但是要加tag
            newForum.insert(items[0]);
        }
    }
    file.close();
    file2.close(); 
    
    // new forum has moderator person
    file2.open(directory1+"/forum_hasModerator_person_0_0.csv",std::ios::app);
    for(auto i=newForum.begin();i!=newForum.end();i++) {
    	std::string forum = *i;
    	for (auto j : forum2person[*i]) {
            if(twoHopResults.find(j)==twoHopResults.end()){
                twoHopResults.insert(j);
                file2<<forum<<"|"<<j<<std::endl;
            }
    		// file2<<forum<<"|"<<j<<std::endl;
		}
	}
    file2.close();

    
	// // forum_hasTag_tag
    // file.open(directory+"/forum_hasTag_tag_0_0.csv");
    // file2.open(directory1+"/forum_hasTag_tag_0_0.csv");
    
    // // std::getline(file, line); // Skip header line
    // while (std::getline(file, line)) {
    //     std::stringstream ss(line);
    //     std::string item;
    //     std::vector<std::string> items;

    //     while (std::getline(ss, item, '|')) {
    //         items.push_back(item);
    //     }

    //     if (items.size() != 2) {
    //         std::cerr << "Invalid line: " << line << std::endl;
    //         continue;
    //     }
        
    //     if(forum.find(items[0])!=forum.end()){
    //         file2<<line<<std::endl;
    //         // forum.insert(items[0]); // 这里加不加都无所谓，提前已经完成了person_forum
    //         // newForum.insert(items[0]);
    //     }
    // }
    // file.close();
    // file2.close(); 

	// forum
    file.open(directory+"/forum_0_0.csv");
    file2.open(directory1+"/forum_0_0.csv");
    
    // std::getline(file, line); // Skip header line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }

        // if (items.size() != 2) {
        //     std::cerr << "Invalid line: " << line << std::endl;
        //     continue;
        // }
        
        if(forum.find(items[0])!=forum.end()){
            file2<<line<<std::endl;
            // forum.insert(items[0]); // 这里加不加都无所谓，提前已经完成了person_forum
            // newForum.insert(items[0]);
        }
    }
    file.close();
    file2.close(); 

	// person
    file.open(directory+"/person_0_0.csv");
    file2.open(directory1+"/person_0_0.csv");
    
    // std::getline(file, line); // Skip header line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }

        // if (items.size() != 2) {
        //     std::cerr << "Invalid line: " << line << std::endl;
        //     continue;
        // }
        
        if(twoHopResults.find(items[0])!=twoHopResults.end()){
            file2<<line<<std::endl;
            // forum.insert(items[0]); // 这里加不加都无所谓，提前已经完成了person_forum
            // newForum.insert(items[0]);
        }
    }
    file.close();
    file2.close(); 

	// comment
    file.open(directory+"/comment_0_0.csv");
    file2.open(directory1+"/comment_0_0.csv");
    
    // std::getline(file, line); // Skip header line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }

        // if (items.size() != 2) {
        //     std::cerr << "Invalid line: " << line << std::endl;
        //     continue;
        // }
        
        if(comment.find(items[0])!=comment.end()){
            file2<<line<<std::endl;
            // forum.insert(items[0]); // 这里加不加都无所谓，提前已经完成了person_forum
            // newForum.insert(items[0]);
        }
    }
    file.close();
    file2.close(); 

	// post
    file.open(directory+"/post_0_0.csv");
    file2.open(directory1+"/post_0_0.csv");
    
    // std::getline(file, line); // Skip header line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }

        // if (items.size() != 2) {
        //     std::cerr << "Invalid line: " << line << std::endl;
        //     continue;
        // }
        
        if(post.find(items[0])!=post.end()){
            file2<<line<<std::endl;
            // forum.insert(items[0]); // 这里加不加都无所谓，提前已经完成了person_forum
            // newForum.insert(items[0]);
        }
    }
    file.close();
    file2.close(); 

    std::cout<<"comment num:"<<comment.size()<<std::endl;
    std::cout<<"post num:"<<post.size()<<std::endl;
    std::cout<<"forum num:"<<forum.size()<<std::endl;
       

}

void loadExknows(std::string filename,
                std::map<std::string, std::map<std::string,std::string>>& knows,
                std::unordered_map<std::string,std::unordered_set<std::string>>& person_hasInterest_tag,
                std::unordered_map<std::string,std::unordered_set<std::string>>& person_isLocatedIn_place,
                std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& person_studyAt_organisation,
                std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& person_workAt_organisation,
                std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& forum_hasMember_person_r,
                std::unordered_map<std::string,std::unordered_set<std::string>>& forum_hasModerator_person_r,
                std::unordered_map<std::string,std::unordered_set<std::string>>& forum_hasModerator_person,
                std::unordered_map<std::string,std::unordered_set<std::string>>& comment_hasCreator_person_r,
                std::unordered_map<std::string,std::unordered_set<std::string>>& comment_hasCreator_person,
                std::unordered_map<std::string,std::unordered_set<std::string>>& post_hasCreator_person_r,
                std::unordered_map<std::string,std::unordered_set<std::string>>& post_hasCreator_person,
                std::unordered_map<std::string,std::unordered_set<std::string>>& post_hasTag_tag,
                std::unordered_map<std::string,std::unordered_set<std::string>>& post_isLocatedIn_place,
                std::unordered_map<std::string,std::unordered_set<std::string>>& comment_hasTag_tag,
                std::unordered_map<std::string,std::unordered_set<std::string>>& comment_isLocatedIn_place,
                std::unordered_map<std::string,std::unordered_set<std::string>>& forum_containerOf_post_r,
                std::unordered_map<std::string, std::unordered_set<std::string>>& reply1,
                std::unordered_map<std::string, std::unordered_set<std::string>>& reply2,
                std::unordered_map<std::string, std::string>& person_0_0,
                std::unordered_map<std::string, std::string>& comment_0_0,
                std::unordered_map<std::string, std::string>& post_0_0,
                std::unordered_map<std::string, std::string>& forum_0_0) {
    // 获取文件夹路径
    // std::string directory = filename.substr(0, filename.find_last_of('/'));
    std::string directory = filename;

    std::cout<<"inputdir:"<<directory<<std::endl;
    
    // person_hasInterest_tag
    
    std::ifstream file(directory+"/person_hasInterest_tag_0_0.csv");
    
    std::string line;

    // std::getline(file, line); // Skip header line
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
        person_hasInterest_tag[items[0]].insert(items[1]);
    }
    file.close();
    
    // person_isLocatedIn_city
    
    file.open(directory+"/person_isLocatedIn_place_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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
        person_isLocatedIn_place[items[0]].insert(items[1]);
    }
    file.close();
    
    // person_studyAt_university
    
    file.open(directory+"/person_studyAt_organisation_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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
        person_studyAt_organisation[items[0]][items[1]] = items[2];
    }
    file.close();
    
    
    // person_workAt_company
    
    file.open(directory+"/person_workAt_organisation_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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
        person_workAt_organisation[items[0]][items[1]] = items[2];
    }
    file.close();
    
    
    // forum_hasMember_person
    
    file.open(directory+"/forum_hasMember_person_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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
        forum_hasMember_person_r[items[1]][items[0]] = items[2];
    }
    file.close();
    
    // forum_hasModerator_person

    file.open(directory+"/forum_hasModerator_person_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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
        forum_hasModerator_person_r[items[1]].insert(items[0]);
        forum_hasModerator_person[items[0]].insert(items[1]);
    }
    file.close();

    // comment hasCreator person
    file.open(directory+"/comment_hasCreator_person_0_0.csv");
    // file2<<"Comment.id|Person.id"<<std::endl;    //为什么没有输出到呢，open会覆盖

    // std::getline(file, line); // Skip header line
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
        comment_hasCreator_person[items[0]].insert(items[1]);
    }
    file.close();

    // post hasCreator person
    file.open(directory+"/post_hasCreator_person_0_0.csv");
    // file2<<"Post.id|Person.id"<<std::endl;

    // std::getline(file, line); // Skip header line
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
        post_hasCreator_person[items[0]].insert(items[1]);
    }
    file.close();

    file.open(directory+"/comment_replyOf_comment_0_0.csv");
    // std::getline(file, line); // Skip header line
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

        reply1[items[0]].insert(items[1]);
    }
    file.close();

    file.open(directory+"/comment_replyOf_post_0_0.csv");
    // std::getline(file, line); // Skip header line
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

        reply2[items[0]].insert(items[1]);
    }
    file.close();
    
    
	// post_hasTag_tag
    file.open(directory+"/post_hasTag_tag_0_0.csv");
    // std::set <std::string> newPerson;
    
    // std::getline(file, line); // Skip header line
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

	// post_isLocatedIn_place
    file.open(directory+"/post_isLocatedIn_place_0_0.csv");
    // std::set <std::string> newPerson;
    
    // std::getline(file, line); // Skip header line
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
    
    // std::getline(file, line); // Skip header line
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
    
    // std::getline(file, line); // Skip header line
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
    
    
	// forum_containerOf_post
    file.open(directory+"/forum_containerOf_post_0_0.csv");
    std::set <std::string> newForum;
    
    // std::getline(file, line); // Skip header line
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

	// forum
    file.open(directory+"/forum_0_0.csv");
    
    // std::getline(file, line); // Skip header line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }
        forum_0_0[items[0]] = line;
    }
    file.close();

	// person
    file.open(directory+"/person_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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

	// comment
    file.open(directory+"/comment_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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
    
    // std::getline(file, line); // Skip header line
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

}

void divideInMemory(std::set<std::string> &twoHopResults, std::string outdir, std::set<std::string>& comment, std::set<std::string>& post, 
                std::map<std::string, std::map<std::string,std::string>>& knows,
                std::unordered_map<std::string,std::unordered_set<std::string>>& person_hasInterest_tag,
                std::unordered_map<std::string,std::unordered_set<std::string>>& person_isLocatedIn_place,
                std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& person_studyAt_organisation,
                std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& person_workAt_organisation,
                std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& forum_hasMember_person_r,
                std::unordered_map<std::string,std::unordered_set<std::string>>& forum_hasModerator_person_r,
                std::unordered_map<std::string,std::unordered_set<std::string>>& forum_hasModerator_person,
                std::unordered_map<std::string,std::unordered_set<std::string>>& comment_hasCreator_person_r,
                std::unordered_map<std::string,std::unordered_set<std::string>>& comment_hasCreator_person,
                std::unordered_map<std::string,std::unordered_set<std::string>>& post_hasCreator_person_r,
                std::unordered_map<std::string,std::unordered_set<std::string>>& post_hasCreator_person,
                std::unordered_map<std::string,std::unordered_set<std::string>>& post_hasTag_tag,
                std::unordered_map<std::string,std::unordered_set<std::string>>& post_isLocatedIn_place,
                std::unordered_map<std::string,std::unordered_set<std::string>>& comment_hasTag_tag,
                std::unordered_map<std::string,std::unordered_set<std::string>>& comment_isLocatedIn_place,
                std::unordered_map<std::string,std::unordered_set<std::string>>& forum_containerOf_post_r,
                std::unordered_map<std::string, std::unordered_set<std::string>>& reply1,
                std::unordered_map<std::string, std::unordered_set<std::string>>& reply2,
                std::unordered_map<std::string, std::string>& person_0_0,
                std::unordered_map<std::string, std::string>& comment_0_0,
                std::unordered_map<std::string, std::string>& post_0_0,
                std::unordered_map<std::string, std::string>& forum_0_0) {
    // 获取文件夹路径
    // std::string directory = filename.substr(0, filename.find_last_of('/'));
    std::string directory1 = "./0/";
    directory1 = outdir;

    // 如果文件夹不存在，则创建它
    if (!createDirectory(directory1)) {
        std::cerr << "Failed or no need to create directory." << std::endl;
    }
    std::cout<<"outputdir:"<<directory1<<std::endl;
    
    std::set<std::string> forum;

    // person_hasInterest_tag
    std::ofstream file2(directory1+"/person_hasInterest_tag_0_0.csv");
    for(auto i:twoHopResults){
        if(person_hasInterest_tag.find(i)!=person_hasInterest_tag.end()){
            for(auto j:person_hasInterest_tag[i]){
                file2<<i<<'|'<<j<<std::endl;
            }
        }
    }
    file2.close(); 
    
    // person_isLocatedIn_city
    
    file2.open(directory1+"/person_isLocatedIn_place_0_0.csv");
    
    for(auto i:twoHopResults){
        if(person_isLocatedIn_place.find(i)!=person_isLocatedIn_place.end()){
            for(auto j:person_isLocatedIn_place[i]){
                file2<<i<<'|'<<j<<std::endl;
            }
        }
    }
    file2.close(); 
    
    // person_studyAt_university
    
    
    file2.open(directory1+"/person_studyAt_organisation_0_0.csv");
    
    for(auto i:twoHopResults){
        if(person_studyAt_organisation.find(i)!=person_studyAt_organisation.end()){
            for(auto j:person_studyAt_organisation[i]){
                file2<<i<<'|'<<j.first<<'|'<<j.second<<std::endl;
            }
        }
    }
    file2.close(); 
    
    // person_workAt_company

    file2.open(directory1+"/person_workAt_organisation_0_0.csv");
    
    for(auto i:twoHopResults){
        if(person_workAt_organisation.find(i)!=person_workAt_organisation.end()){
            for(auto j:person_workAt_organisation[i]){
                file2<<i<<'|'<<j.first<<'|'<<j.second<<std::endl;
            }
        }
    }
    file2.close(); 
    
    
    // forum_hasMember_person
    
    file2.open(directory1+"/forum_hasMember_person_0_0.csv");
    for(auto i:twoHopResults){
        if(forum_hasMember_person_r.find(i)!=forum_hasMember_person_r.end()){
            for(auto j:forum_hasMember_person_r[i]){
                file2<<j.first<<'|'<<i<<'|'<<j.second<<std::endl;
                forum.insert(j.first);
            }
        }
    }
    file2.close(); 
    
    // forum_hasModerator_person

    file2.open(directory1+"/forum_hasModerator_person_0_0.csv");
    for(auto i:twoHopResults){
        if(forum_hasModerator_person_r.find(i)!=forum_hasModerator_person_r.end()){
            for(auto j:forum_hasModerator_person_r[i]){
                file2<<j<<'|'<<i<<std::endl;
                forum.insert(j);
            }
        }
    }
    file2.close(); 

    // comment hasCreator person
    file2.open(directory1+"/comment_hasCreator_person_0_0.csv");
    for(auto i:twoHopResults){
        if(comment_hasCreator_person_r.find(i)!=comment_hasCreator_person_r.end()){
            for(auto j:comment_hasCreator_person_r[i]){
                // file2<<j<<'|'<<i<<std::endl;
                comment.insert(j);
            }
        }
    }
    file2.close();

    // post hasCreator person
    file2.open(directory1+"/post_hasCreator_person_0_0.csv");
    
    for(auto i:twoHopResults){
        if(post_hasCreator_person_r.find(i)!=post_hasCreator_person_r.end()){
            for(auto j:post_hasCreator_person_r[i]){
                file2<<j<<'|'<<i<<std::endl;
                post.insert(j);
            }
        }
    }
    file2.close();
    

	// commet链至post
    // dfs+路径压缩
    std::ofstream file3_1(directory1+"/comment_replyOf_comment_0_0.csv");
    std::ofstream file3_2(directory1+"/comment_replyOf_post_0_0.csv");
    // 去重comment2comment
    std::unordered_map<std::string,std::unordered_map<std::string,int>> visited; 
    // 去重comment2post
    std::unordered_map<std::string,std::unordered_map<std::string,int>> visited1; 

    // dfs得到 所有的commentPath(comment commment-neigbor comment-end)
    std::unordered_map<std::string,std::vector<std::string>> commentPath;
    for (auto i=comment.begin();i!=comment.end();i++) {
        std::vector<std::string> tmp;
        dfs(*i, reply1, tmp, commentPath);
	}
    // 路径压缩
    // 从引入的comment开始dfs得到（起点comment，一跳comment，终点comment）
    // 1.（起点，，）
    // 加入起点comment的post的边及邻居
    // 2.（起点，一跳，）
    // 加入起点到一跳的边及邻居，加入起点到【一跳comment的post邻居】和起点到【一跳comment的post邻居】边
    // 2.（起点，一跳，终点）
    // 加入起点到一跳的边及邻居，加入起点到【终点comment的post邻居】和起点到【终点comment的post邻居】边
    for(auto i=commentPath.begin();i!=commentPath.end();i++) {
        std::string from = i->first;
    	std::vector<std::string> tmp = i->second;
        if(tmp.size()==1) {
            // comment没有comment邻居
            // ！！！这里也要存post
            comment.insert(tmp[0]); // 里面本来就存了
            // file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
            for(std::string j : reply2[tmp[0]]) {
                post.insert(j);
                if(visited1[tmp[0]].find(j)==visited1[tmp[0]].end()) {
                    visited1[tmp[0]][j] = 1;
                    file3_2<<tmp[0]<<'|'<<j<<std::endl;
                }
                // file3_2<<tmp[0]<<'|'<<j<<std::endl;
            }
        } else if(tmp.size()==2) {
            // comment只有一步comment
            comment.insert(tmp[1]);
            if(visited[tmp[0]].find(tmp[1])==visited[tmp[0]].end()) {
                visited[tmp[0]][tmp[1]] = 1;
                file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
            }
            // file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
            for(std::string j : reply2[tmp[1]]) {
                // comment存在comment链
                post.insert(j);
                if(visited1[tmp[1]].find(j)==visited1[tmp[1]].end()) {
                    visited1[tmp[1]][j] = 1;
                    file3_2<<tmp[1]<<'|'<<j<<std::endl;
                }
                // file3_2<<tmp[1]<<'|'<<j<<std::endl;
            }
        } else {
            // 顶点有多步邻居
            comment.insert(tmp[1]);
            if(visited[tmp[0]].find(tmp[1])==visited[tmp[0]].end()) {
                visited[tmp[0]][tmp[1]] = 1;
                file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
            }
            // file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
            // 这里把末端邻居合并为post
            for(std::string j : reply2[tmp[2]]) {
                post.insert(j);
                if(visited1[tmp[1]].find(j)==visited1[tmp[1]].end()) {
                    visited1[tmp[1]][j] = 1;
                    file3_2<<tmp[1]<<'|'<<j<<std::endl;
                }
                // 压缩了路径，不过好像导致comment2post变多了？
            }
        }
    }
	
	// comment_hasCreator_person
    file2.open(directory1+"/comment_hasCreator_person_0_0.csv");// 注意这里防止覆盖app打开
    // 这里有所有的新的comment，这里所有的comment到person的边都是需要加的，所以前面不需要加，这里就是要覆盖
    // 新comment到person
    std::set <std::string> newPerson;
    for(auto c:comment) {
        if(comment_hasCreator_person.find(c)!=comment_hasCreator_person.end()) {
            for(auto p:comment_hasCreator_person[c]) {
                if(twoHopResults.find(p)==twoHopResults.end()) {
                    newPerson.insert(p);    //newPerson干啥了？好像没用吧
                    twoHopResults.insert(p);
                }
                file2<<c<<'|'<<p<<std::endl;    // 所有comment到person的边都是需要加的
            }
        }
    }
    file2.close();
    
    
	// post_hasCreator_person
    file2.open(directory1+"/post_hasCreator_person_0_0.csv", std::ios::app);
    for(auto p:post) {
        if(post_hasCreator_person.find(p)!=post_hasCreator_person.end()) {
            for(auto c:post_hasCreator_person[p]) {
                if(twoHopResults.find(c)==twoHopResults.end()) {
                    newPerson.insert(c);
                    twoHopResults.insert(c);
                    file2<<p<<'|'<<c<<std::endl;    // 不是新person不用再引入
                }
            }
        }
    }
    file2.close();
    
    // todo:newPerson与其他person的knows,finished!
    // 是否需要传入knows? 传入吧 直接打开新建也是一样
    // 遍历twoHopResults中的每个person
    file2.open(directory1+"/person_knows_person_0_0.csv");
    std::map<std::string, std::map<std::string,int>> kvisited;
    for (const auto& person : twoHopResults) {
        if (knows.find(person) != knows.end()) {
            // 找到person的所有邻居
            for (const auto& neighbor : knows[person]) {
                // 将每个person与其邻居之间的knows边存储起来
                if(twoHopResults.find(neighbor.first)!=twoHopResults.end()){
                    if(!isEdgeVisited(neighbor.first, person, kvisited)) {
                        kvisited[person][neighbor.first] = 1;
                        file2<<person<<'|'<<neighbor.first<<'|'<<neighbor.second<<std::endl;
                    }
                }
            }
        }
    }
    file2.close();
    
    
	// post_hasTag_tag
    file2.open(directory1+"/post_hasTag_tag_0_0.csv");
    for(auto p:post) {
        if(post_hasTag_tag.find(p)!=post_hasTag_tag.end()) {
            for(auto t:post_hasTag_tag[p]) {
                file2<<p<<'|'<<t<<std::endl;
            }
        }
    }
    file2.close(); 

	// post_isLocatedIn_place
    file2.open(directory1+"/post_isLocatedIn_place_0_0.csv");
    for(auto p:post) {
        if(post_isLocatedIn_place.find(p)!=post_isLocatedIn_place.end()) {
            for(auto t:post_isLocatedIn_place[p]) {
                file2<<p<<'|'<<t<<std::endl;
            }
        }
    }
    file2.close(); 
	
	
	// comment_hasTag_tag
    file2.open(directory1+"/comment_hasTag_tag_0_0.csv");
    for(auto c:comment) {
        if(comment_hasTag_tag.find(c)!=comment_hasTag_tag.end()) {
            for(auto t:comment_hasTag_tag[c]) {
                file2<<c<<'|'<<t<<std::endl;
            }
        }
    }
    file2.close(); 
    
	// comment_isLocatedIn_place
    file2.open(directory1+"/comment_isLocatedIn_place_0_0.csv");
    for(auto c:comment) {
        if(comment_isLocatedIn_place.find(c)!=comment_isLocatedIn_place.end()) {
            for(auto t:comment_isLocatedIn_place[c]) {
                file2<<c<<'|'<<t<<std::endl;
            }
        }
    }
    file2.close(); 
    
    
	// forum_containerOf_post
    file2.open(directory1+"/forum_containerOf_post_0_0.csv");
    std::set <std::string> newForum;
    for(auto p: post) {
        if(forum_containerOf_post_r.find(p)!=forum_containerOf_post_r.end()) {
            for(auto f:forum_containerOf_post_r[p]) {
                file2<<f<<'|'<<p<<std::endl;
                forum.insert(f);
                if(forum.find(f)==forum.end()) {
                    newForum.insert(f);
                }
            }
        }
    }
    file2.close(); 
    
    // new forum has moderator person
    file2.open(directory1+"/forum_hasModerator_person_0_0.csv",std::ios::app);
    for(auto i=newForum.begin();i!=newForum.end();i++) {
    	std::string forum = *i;
    	for (auto j : forum_hasModerator_person[forum]) {
            if(twoHopResults.find(j)==twoHopResults.end()){
                twoHopResults.insert(j);
                file2<<forum<<"|"<<j<<std::endl;    //已经存在的人是不会需要加这条边的
            }
		}
	}
    file2.close();

    // 不需要这个！
    // file2.open(directory1+"/forum_hasTag_tag_0_0.csv");
    // for(auto i:forum){
    //     if(forum_hasTag_tag.find(i)!=forum_hasTag_tag.end()){
    //         for(auto j:forum_hasTag_tag[i]){
    //             file2<<i<<'|'<<j<<std::endl;
    //         }
    //     }
    // }
    // file2.close(); 

	// forum
    file2.open(directory1+"/forum_0_0.csv");
    for(auto i:forum){
        file2<<forum_0_0[i]<<std::endl;
    }
    file2.close(); 

	// person
    file2.open(directory1+"/person_0_0.csv");
    for(auto i:twoHopResults){
        file2<<person_0_0[i]<<std::endl;
    }
    file2.close(); 

	// comment
    file2.open(directory1+"/comment_0_0.csv");
    for(auto i:comment){
        file2<<comment_0_0[i]<<std::endl;
    }
    file2.close(); 

	// post
    file2.open(directory1+"/post_0_0.csv");
    for(auto i:post){
        file2<<post_0_0[i]<<std::endl;
    }
    file2.close(); 

    std::cout<<"comment num:"<<comment.size()<<std::endl;
    std::cout<<"post num:"<<post.size()<<std::endl;
    std::cout<<"forum num:"<<forum.size()<<std::endl;
       

}

// 函数用于计算每个 person 的邻居数目并存储在 neighbor_counts 中
void calculateNeighborCounts(const std::map<std::string, Person>& people,
                              const std::map<std::string, std::map<std::string, std::string>>& knows,
                              std::map<std::string, int>& neighbor_counts) {
    // 初始化邻居计数映射，为每个 person 设置计数为 0
    neighbor_counts.clear();

    // 遍历 people 映射中的所有 person
    for (const auto& person_pair : people) {
        const std::string& person_id = person_pair.first;
        // 为每个人初始化邻居计数
        neighbor_counts[person_id] = 0;
    }

    // 再次遍历 people 映射中的所有 person
    // 对于每个人，检查 knows 映射中是否存在他们的邻居
    for (const auto& person_pair : people) {
        const std::string& person_id = person_pair.first;
        auto knows_it = knows.find(person_id);
        if (knows_it != knows.end()) {
            // 如果 knows 中有这个 person 的记录，统计他们的邻居数目
            neighbor_counts[person_id] = knows_it->second.size();
        }
    }
}


std::map<std::string, int> calculateTwoHopNeighborCounts(
    const std::map<std::string, Person>& people,
    const std::map<std::string, std::map<std::string, std::string>>& knows,
    std::map<std::string, std::set<std::string>>& two_hop_people) {
    
    std::map<std::string, int> two_hop_neighbor_counts;

    // 初始化two_hop_neighbor_counts为所有people中的person的两跳邻居数目为0
    for (const auto& person_pair : people) {
        const std::string& person_id = person_pair.first;
        two_hop_neighbor_counts[person_id] = 0;
        two_hop_people[person_id]; // 初始化集合
    }

    // 遍历knows映射，统计两跳邻居数目
    for (const auto& outer_pair : knows) {
        const std::string& person = outer_pair.first;
        for (const auto& inner_pair : outer_pair.second) {
            const std::string& one_hop_neighbor = inner_pair.first;
            two_hop_people[person].insert(one_hop_neighbor);
            // 检查一跳邻居是否有两跳邻居
            auto two_hop_it = knows.find(one_hop_neighbor);
            if (two_hop_it != knows.end()) {
                // 更新两跳邻居集合
                for (const auto& two_hop_pair : two_hop_it->second) {
                    const std::string& two_hop_neighbor = two_hop_pair.first;
                    two_hop_people[person].insert(two_hop_neighbor);
                }
            }
        }
        // 更新两跳邻居数目
        two_hop_neighbor_counts[person] += two_hop_people[person].size();
    }

    return two_hop_neighbor_counts;
}

// int calculateSetTwoHopNeighborNum(const std::set<std::string>& people, const std::map<std::string, int>& two_hop_people) {
//     int num = 0;
//     for (const std::string& person : people) {
//         auto it = two_hop_people.find(person);
//         if (it != two_hop_people.end()) { 
//             num += it->second;
//         }
//     }
//     return num;
// }

int calculateSetTwoHopNeighborNum(const std::set<std::string>& center, const std::map<std::string, int>& two_hop_neighbor_counts) {
    int total = 0;
    for (const auto& person : center) {
        auto it = two_hop_neighbor_counts.find(person);
        if (it != two_hop_neighbor_counts.end()) {
            total += it->second;
        } else {
            std::cerr << "Person " << person << " not found in two_hop_neighbor_counts" << std::endl;
        }
    }
    return total;
}

int calculateSetTwoHopNeighborNumNew(int oldthreshhold, std::string neighbor, std::set<std::string> twoHop, std::map<std::string, std::set<std::string>> two_hop_neighbors) {
    // 这里不能把新的邻居加入到twoHop中，因为阈值不超过就不能使用新的两跳邻居列表
    int total = oldthreshhold;
    for (const auto& person : two_hop_neighbors[neighbor]) {
        if(twoHop.find(person)==twoHop.end())
            total++;
    }
    return total;
}

std::set<std::string> getKnowsNeighbor(const std::set<std::string>& persons, 
                                        const std::map<std::string, std::map<std::string, std::string>>& knows) {
    std::set<std::string> neighbor;
    
    // 遍历传入的persons集合
    for (const std::string& person : persons) {
        // 在knows映射中查找当前person的邻居
        auto knows_it = knows.find(person);
        if (knows_it != knows.end()) {
            // 遍历并添加每个邻居到neighbor集合中
            for (const auto& neighbor_pair : knows_it->second) {
                neighbor.insert(neighbor_pair.first); // 假设neighbor_pair.first是邻居的ID
            }
        }
    }
    
    return neighbor;
}

std::set<std::string> getKnowsNeighbor1(const std::set<std::string>& persons, 
                                        const std::map<std::string, std::map<std::string, std::string>>& knows, const std::set<std::string>& peopleSet) {
    std::set<std::string> neighbor;
    
    // 遍历传入的persons集合
    for (const std::string& person : persons) {
        // 在knows映射中查找当前person的邻居
        auto knows_it = knows.find(person);
        if (knows_it != knows.end()) {
            // 遍历并添加每个邻居到neighbor集合中
            for (const auto& neighbor_pair : knows_it->second) {
                // 过滤已经作为center的person
                if(peopleSet.find(neighbor_pair.first)!=peopleSet.end())
                    neighbor.insert(neighbor_pair.first); // 假设neighbor_pair.first是邻居的ID
            }
        }
    }
    
    return neighbor;
}

std::string getSmallestOneHop(std::set<std::string>& persons, std::map<std::string, int>& hop) {
    if (persons.empty()) {
        return ""; // 如果persons集合为空，返回空字符串
    }

    // 初始化最小一跳邻居数目为最大整数值
    int smallest = 99999999;
    std::string personWithFewestNeighbors;

    // 遍历persons集合，找出具有最小一跳/二跳邻居数目的person
    for (const std::string& person : persons) {
        auto oneHopIt = hop.find(person);
        if (oneHopIt != hop.end()) {
            if (oneHopIt->second < smallest) {
                smallest = oneHopIt->second;
                personWithFewestNeighbors = person;
                // std::cout<<smallest<<std::endl;
                // std::cout<<"small:"<<personWithFewestNeighbors<<std::endl;
            }
        }
    }
    // std::cout<<smallest<<std::endl;
    return personWithFewestNeighbors;
}


std::string getBiggestOneHop(std::set<std::string>& persons, std::map<std::string, int>& hop) {
    if (persons.empty()) {
        return ""; // 如果persons集合为空，返回空字符串
    }

    // 初始化最小一跳邻居数目为最大整数值
    int biggest = -1;
    std::string personWithFewestNeighbors;

    // 遍历persons集合，找出具有最小一跳/二跳邻居数目的person
    for (const std::string& person : persons) {
        auto oneHopIt = hop.find(person);
        if (oneHopIt != hop.end()) {
            if (oneHopIt->second > biggest) {
                biggest = oneHopIt->second;
                personWithFewestNeighbors = person;
                // std::cout<<biggest<<std::endl;
                // std::cout<<"small:"<<personWithFewestNeighbors<<std::endl;
            }
        }
    }
    // std::cout<<biggest<<std::endl;
    return personWithFewestNeighbors;
}

// 计算两个集合的 Jaccard 相似度
double jaccard_similarity(const std::set<std::string>& set1, const std::set<std::string>& set2) {
    std::set<std::string> intersection;
    std::set<std::string> union_set;

    // 计算交集
    std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), std::inserter(intersection, intersection.begin()));
    // 计算并集
    std::set_union(set1.begin(), set1.end(), set2.begin(), set2.end(), std::inserter(union_set, union_set.begin()));

    return static_cast<double>(intersection.size()) / union_set.size();
}

// 合并两个集合
std::set<std::string> merge_sets(const std::set<std::string>& set1, const std::set<std::string>& set2) {
    std::set<std::string> merged_set = set1;
    merged_set.insert(set2.begin(), set2.end());
    return merged_set;
}

// 合并相似度高的集合
std::set<std::set<std::string>> merge_high_similarity_sets(std::set<std::set<std::string>>& twoHopSetsInput, double threshold) {
    // Step 1: 合并所有孤立的集合
    std::set<std::string> single_value_set;
    std::set<std::set<std::string>> remaining_sets;

    for (const auto& s : twoHopSetsInput) {
        if (s.size() == 1) {
            single_value_set.insert(*s.begin());
        } else {
            remaining_sets.insert(s);
        }
    }

    if (!single_value_set.empty()) {
        remaining_sets.insert(single_value_set);
    }

    // Step 2: 按相似度合并集合
    std::set<std::set<std::string>> twoHopSets = remaining_sets;
    bool merged = true;

    while (merged) {
        merged = false;
        std::set<std::set<std::string>> new_twoHopSets;
        auto it1 = twoHopSets.begin();

        while (it1 != twoHopSets.end()) {
            bool is_merged = false;
            auto it2 = std::next(it1);

            while (it2 != twoHopSets.end()) {
                if (jaccard_similarity(*it1, *it2) >= threshold) {
                    std::set<std::string> mset = merge_sets(*it1, *it2);
                    new_twoHopSets.insert(mset);

                    // 合并后删除原集合
                    it2 = twoHopSets.erase(it2);
                    it1 = twoHopSets.erase(it1);
                    merged = true;
                    is_merged = true;
                    break;
                } else {
                    ++it2;
                }
            }

            if (!is_merged) {
                new_twoHopSets.insert(*it1);
                ++it1;
            }
        }
        twoHopSets = new_twoHopSets;
    }

    return twoHopSets;
}

void partition1(std::string inputdir, std::map<std::string, Person> people, std::map<std::string, std::map<std::string,std::string>> knows,std::map<std::string, std::set<std::string>> two_hop_neighbors) {
    
    std::ofstream outfile("./route/message2person.txt");
    std::ifstream infile(inputdir+"comment_hasCreator_person.csv");
    std::string line;
    while (std::getline(infile, line)) {
        outfile << line << std::endl;
    }
    infile.close();
    infile.open(inputdir+"post_hasCreator_person.csv");
    std::string tmplines;
    while (std::getline(infile, tmplines)) {
        outfile << line << std::endl;
    }
    outfile.close();

    std::ofstream outfile1("./route/person2part.txt");
    for(auto i:people) {
        outfile1<<i.first<<'|'<<i.first<<std::endl;
    }
    outfile1.close();

    
    // 划分方法1：每个person做中心person
    int part = 0;
    for(auto it=people.begin();it!=people.end();it++){
        std::cout<<part<<std::endl;
        
        std::string givenPersonID = it->first; // Example ID
        std::set<std::string> twoHopResults;

        twoHopResults = two_hop_neighbors[givenPersonID];
        twoHopResults.insert(givenPersonID);

        std::cout<<"center PersonID:"<<givenPersonID<<std::endl;

        std::cout<<"twohop persons num:"<<twoHopResults.size()<<std::endl;
        // 邻居最小一跳/二跳，基本上是一致的
        std::set<std::string> center;
        center.insert(givenPersonID);
        std::set<std::string> cNeighbor = getKnowsNeighbor(center,knows);

        std::set<std::string> comment;
        std::set<std::string> post;

        getDivide(twoHopResults, inputdir, "./output/"+givenPersonID+"/", comment, post, knows);
        
        std::cout<<"part persons num:"<<twoHopResults.size()<<std::endl<<std::endl;
        part++;
    }
}

void partition1_1(std::string inputdir, std::map<std::string, Person> people, std::map<std::string, std::map<std::string,std::string>> knows,std::map<std::string, std::set<std::string>> two_hop_neighbors) {
    // 划分方法1：一次load
    std::unordered_map<std::string,std::unordered_set<std::string>> person_hasInterest_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> person_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_studyAt_organisation;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_workAt_organisation;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> forum_hasMember_person_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> forum_hasModerator_person_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> forum_hasModerator_person;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasCreator_person_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasCreator_person;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasCreator_person_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasCreator_person;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasTag_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasTag_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_set<std::string>> forum_containerOf_post_r;
    std::unordered_map<std::string, std::unordered_set<std::string>> reply1;
    std::unordered_map<std::string, std::unordered_set<std::string>> reply2;
    std::unordered_map<std::string, std::string> person_0_0;
    std::unordered_map<std::string, std::string> comment_0_0;
    std::unordered_map<std::string, std::string> post_0_0;
    std::unordered_map<std::string, std::string> forum_0_0;

    loadExknows(inputdir, knows, person_hasInterest_tag,person_isLocatedIn_place,person_studyAt_organisation,person_workAt_organisation,forum_hasMember_person_r,forum_hasModerator_person_r,forum_hasModerator_person,comment_hasCreator_person_r,comment_hasCreator_person,post_hasCreator_person_r,post_hasCreator_person,post_hasTag_tag,post_isLocatedIn_place,comment_hasTag_tag,comment_isLocatedIn_place,forum_containerOf_post_r,reply1,reply2,person_0_0,comment_0_0,post_0_0,forum_0_0);
    
    
    std::ofstream outfile("./route/message2person.txt");
    for(auto i:comment_hasCreator_person){
        for(auto j:comment_hasCreator_person[i.first]){
            outfile<<i.first<<'|'<<j<<std::endl;
        }
    }
    for(auto i:post_hasCreator_person){
        for(auto j:post_hasCreator_person[i.first]){
            outfile<<i.first<<'|'<<j<<std::endl;
        }
    }
    outfile.close();

    std::ofstream outfile1("./route/person2part.txt");
    for(auto i:people) {
        outfile1<<i.first<<'|'<<i.first<<std::endl;
    }
    outfile1.close();

    // 划分方法1.1：每个person做中心person
    int part = 0;
    
    for(auto it=people.begin();it!=people.end();it++){
        std::cout<<part<<std::endl;
        
        std::string givenPersonID = it->first; // Example ID
        std::set<std::string> twoHopResults;

        twoHopResults = two_hop_neighbors[givenPersonID];
        twoHopResults.insert(givenPersonID);

        std::cout<<"center PersonID:"<<givenPersonID<<std::endl;

        std::cout<<"twohop persons num:"<<twoHopResults.size()<<std::endl;
        // 邻居最小一跳/二跳，基本上是一致的
        std::set<std::string> center;
        center.insert(givenPersonID);
        std::set<std::string> cNeighbor = getKnowsNeighbor(center,knows);

        std::set<std::string> comment;
        std::set<std::string> post;

        divideInMemory(twoHopResults, "./output/"+givenPersonID+"/", comment, post, knows, person_hasInterest_tag, person_isLocatedIn_place, person_studyAt_organisation, person_workAt_organisation, forum_hasMember_person_r, forum_hasModerator_person_r, forum_hasModerator_person, comment_hasCreator_person_r, comment_hasCreator_person, post_hasCreator_person_r, post_hasCreator_person, post_hasTag_tag, post_isLocatedIn_place, comment_hasTag_tag, comment_isLocatedIn_place, forum_containerOf_post_r, reply1, reply2, person_0_0, comment_0_0, post_0_0, forum_0_0);
        std::cout<<"part persons num:"<<twoHopResults.size()<<std::endl<<std::endl;
        part++;
    }
}

void partition2(int input_threshhold, std::string inputdir, std::map<std::string, Person> people, std::map<std::string, std::map<std::string,std::string>> knows,std::map<std::string, std::set<std::string>> two_hop_neighbors, std::map<std::string, int> neighbor_counts) {
    // 划分方法2：贪心算法合并person
    int part = 0;
    std::set<std::string> peopleSet;
    for(auto it=people.begin();it!=people.end();it++){
        peopleSet.insert(it->first);
    }
    std::ofstream tmpfile("./route/person2part.txt");
    
    while(peopleSet.size()!=0) {
        std::cout<<"part:"<<part<<std::endl;
        std::string seed = getBiggestOneHop(peopleSet,neighbor_counts);
        if(seed=="") continue;
        
        std::set<std::string> twoHopResults;
        twoHopResults.insert(seed);
        twoHopResults.insert(two_hop_neighbors[seed].begin(), two_hop_neighbors[seed].end());
        
        std::set<std::string> center;
        center.insert(seed);
        std::set<std::string> cNeighbor = getKnowsNeighbor1(center,knows,peopleSet);
        
        int threshhold = twoHopResults.size();
        for (const auto& person : center) {
            cNeighbor.erase(person);
        }
        while(cNeighbor.size()!=0&&threshhold<input_threshhold) {
            
            std::string smallestOneHop = getSmallestOneHop(cNeighbor,neighbor_counts);
            if(smallestOneHop=="") break;
            threshhold = calculateSetTwoHopNeighborNumNew(threshhold, smallestOneHop, twoHopResults, two_hop_neighbors); 
            if(threshhold>1355) break;
            
            center.insert(smallestOneHop);
            twoHopResults.insert(smallestOneHop);
            twoHopResults.insert(two_hop_neighbors[smallestOneHop].begin(), two_hop_neighbors[smallestOneHop].end());
            
           
            // 获取center的邻居
            cNeighbor = getKnowsNeighbor1(center,knows,peopleSet);
            // std::cout<<"cNeighborsize:"<<cNeighbor.size()<<std::endl;
            for (const auto& person : center) {
                cNeighbor.erase(person);
            }
            
        }
        std::cout<<"center persons num:"<<center.size()<<std::endl;
        std::cout<<"threshhold:"<<threshhold<<std::endl;
        std::cout<<"twohop persons num:"<<twoHopResults.size()<<std::endl;
        std::string outPath = "./output/";
        outPath+=std::to_string(part);
        for(auto it=center.begin();it!=center.end();it++){
            tmpfile<<*it<<'|'<<part<<std::endl;
        }

        std::set<std::string> comment;
        std::set<std::string> post;
        getDivide(twoHopResults, inputdir, outPath, comment, post, knows);
        
        std::cout<<"part persons num:"<<twoHopResults.size()<<std::endl<<std::endl;
        
        part++;

        for (const auto& person : center) {
            peopleSet.erase(person);
        }
    }
}

void partition2_1(int input_threshhold, std::string inputdir, std::map<std::string, Person> people, std::map<std::string, std::map<std::string,std::string>> knows,std::map<std::string, std::set<std::string>> two_hop_neighbors, std::map<std::string, int> neighbor_counts) {
     // 划分方法2.1：一次load 贪心算法合并person
    std::unordered_map<std::string,std::unordered_set<std::string>> person_hasInterest_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> person_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_studyAt_organisation;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_workAt_organisation;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> forum_hasMember_person_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> forum_hasModerator_person_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> forum_hasModerator_person;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasCreator_person_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasCreator_person;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasCreator_person_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasCreator_person;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasTag_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasTag_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_set<std::string>> forum_containerOf_post_r;
    std::unordered_map<std::string, std::unordered_set<std::string>> reply1;
    std::unordered_map<std::string, std::unordered_set<std::string>> reply2;
    std::unordered_map<std::string, std::string> person_0_0;
    std::unordered_map<std::string, std::string> comment_0_0;
    std::unordered_map<std::string, std::string> post_0_0;
    std::unordered_map<std::string, std::string> forum_0_0;

    loadExknows(inputdir,knows, person_hasInterest_tag,person_isLocatedIn_place,person_studyAt_organisation,person_workAt_organisation,forum_hasMember_person_r,forum_hasModerator_person_r,forum_hasModerator_person,comment_hasCreator_person_r,comment_hasCreator_person,post_hasCreator_person_r,post_hasCreator_person,post_hasTag_tag,post_isLocatedIn_place,comment_hasTag_tag,comment_isLocatedIn_place,forum_containerOf_post_r,reply1,reply2,person_0_0,comment_0_0,post_0_0,forum_0_0);
    
    std::ofstream outfile("./route/message2person.txt");
    for(auto i:comment_hasCreator_person){
        for(auto j:comment_hasCreator_person[i.first]){
            outfile<<i.first<<'|'<<j<<std::endl;
        }
    }
    for(auto i:post_hasCreator_person){
        for(auto j:post_hasCreator_person[i.first]){
            outfile<<i.first<<'|'<<j<<std::endl;
        }
    }
    outfile.close();

    int part = 0;
    std::set<std::string> peopleSet;
    for(auto it=people.begin();it!=people.end();it++){
        peopleSet.insert(it->first);
    }
    std::ofstream tmpfile("./route/person2part.txt");
    
    while(peopleSet.size()!=0) {
        std::cout<<"part:"<<part<<std::endl;
        std::string seed = getBiggestOneHop(peopleSet,neighbor_counts);
        if(seed=="") continue;
        
        std::set<std::string> twoHopResults;
        twoHopResults.insert(seed);
        twoHopResults.insert(two_hop_neighbors[seed].begin(), two_hop_neighbors[seed].end());
        
        std::set<std::string> center;
        center.insert(seed);
        std::set<std::string> cNeighbor = getKnowsNeighbor1(center,knows,peopleSet);
        
        int threshhold = twoHopResults.size();
        for (const auto& person : center) {
            cNeighbor.erase(person);
        }
        while(cNeighbor.size()!=0&&threshhold<input_threshhold) {
            
            std::string smallestOneHop = getSmallestOneHop(cNeighbor,neighbor_counts);
            if(smallestOneHop=="") break;
            threshhold = calculateSetTwoHopNeighborNumNew(threshhold, smallestOneHop, twoHopResults, two_hop_neighbors); 
            // if(threshhold>1355) break;
            
            center.insert(smallestOneHop);
            twoHopResults.insert(smallestOneHop);
            twoHopResults.insert(two_hop_neighbors[smallestOneHop].begin(), two_hop_neighbors[smallestOneHop].end());
            
           
            // 获取center的邻居
            cNeighbor = getKnowsNeighbor1(center,knows,peopleSet);
            // std::cout<<"cNeighborsize:"<<cNeighbor.size()<<std::endl;
            for (const auto& person : center) {
                cNeighbor.erase(person);
            }
            
        }
        std::cout<<"center persons num:"<<center.size()<<std::endl;
        std::cout<<"threshhold:"<<threshhold<<std::endl;
        std::cout<<"twohop persons num:"<<twoHopResults.size()<<std::endl;
        std::string outPath = "./output/";
        outPath+=std::to_string(part);
        for(auto it=center.begin();it!=center.end();it++){
            tmpfile<<*it<<' '<<part<<std::endl;
        }

        std::set<std::string> comment;
        std::set<std::string> post;
        divideInMemory(twoHopResults, outPath, comment, post, knows, person_hasInterest_tag, person_isLocatedIn_place, person_studyAt_organisation, person_workAt_organisation, forum_hasMember_person_r, forum_hasModerator_person_r, forum_hasModerator_person, comment_hasCreator_person_r, comment_hasCreator_person, post_hasCreator_person_r, post_hasCreator_person, post_hasTag_tag, post_isLocatedIn_place, comment_hasTag_tag, comment_isLocatedIn_place, forum_containerOf_post_r, reply1, reply2, person_0_0, comment_0_0, post_0_0, forum_0_0);
        
        std::cout<<"part persons num:"<<twoHopResults.size()<<std::endl<<std::endl;
        
        part++;

        for (const auto& person : center) {
            peopleSet.erase(person);
        }
    }

}

void partition3(int input_threshhold, int input_similarity, std::string inputdir, std::map<std::string, Person> people, std::map<std::string, std::map<std::string,std::string>> knows,std::map<std::string, std::set<std::string>> two_hop_neighbors, std::map<std::string, int> neighbor_counts) {
    // person2part还有问题！
    // 划分方法3：3.1贪心算法合并person，3.2合并相似度高的两跳person集
    std::set<std::set<std::string>> twoHopSets;   // 记录所有分片的两步邻居

    std::set<std::string> peopleSet;    // 记录剩余未做种子分片的person，初始为所有person
    for(auto it=people.begin();it!=people.end();it++){
        peopleSet.insert(it->first);
    }

    int testCenter=0;
    int testtwohop=0;

    while(peopleSet.size()!=0) {    // 存在未作种子的person
        // 取seed，即一跳邻居最多的person
        std::string seed = getBiggestOneHop(peopleSet,neighbor_counts);
        if(seed=="") continue;
        
        //将seed及其两跳邻居加入分片
        std::set<std::string> twoHopResults;
        twoHopResults.insert(seed);
        twoHopResults.insert(two_hop_neighbors[seed].begin(), two_hop_neighbors[seed].end());
        

        // 将seed加入种子集即center
        std::set<std::string> center;
        center.insert(seed);
        // 获取seed的一跳邻居，即候选种子，不包含已经作为center的person
        std::set<std::string> cNeighbor = getKnowsNeighbor1(center,knows,peopleSet);
        for (const auto& person : center) {
            cNeighbor.erase(person);
        }

        int threshhold = twoHopResults.size();  // 阈值初始为seed的两跳邻居数目
        while(cNeighbor.size()!=0&&threshhold<input_threshhold) {

            // 在邻居中找到一跳邻居最小的person
            std::string smallestOneHop = getSmallestOneHop(cNeighbor,neighbor_counts);
            if (smallestOneHop == "") { // cNeighbor为空时返回为""，正常情况下循环直接结束其实不会执行
                break;
            }
            // 更新阈值
            threshhold = calculateSetTwoHopNeighborNumNew(threshhold, smallestOneHop, twoHopResults, two_hop_neighbors);
            // if(threshhold>1700) break;

            // 符合阈值则加入center和两步邻居
            center.insert(smallestOneHop);
            twoHopResults.insert(smallestOneHop);
            twoHopResults.insert(two_hop_neighbors[smallestOneHop].begin(), two_hop_neighbors[smallestOneHop].end());
            // peopleSet.erase(smallestOneHop); //这里可删可不删，cNeighbor.erase(center);会删掉

            // 更新候选种子集，这里已经在两步邻居里的一步邻居还要做种子，相当于从比较小的扩张
            cNeighbor = getKnowsNeighbor1(center,knows,peopleSet);
            // 保证做过center的不会再被选为候选种子
            for (const auto& person : center) {
                cNeighbor.erase(person);
            }
            
        }
        
        // 插入到twoHopSets中待合并
        twoHopSets.insert(twoHopResults);
        
        // 保证做过种子的不会再在其他分片做种子
        for (const auto& person : center) {
            peopleSet.erase(person);
        }

        testCenter+=center.size();
        testtwohop+=twoHopResults.size();
    }
    std::cout<<"all center person num:"<<testCenter<<std::endl;
    std::cout<<"all twohop person num:"<<testtwohop<<std::endl;
    
    // 合并相似度高的两跳person集
    double merge_threshold = input_similarity;  // 相似度阈值,为0时合为一个分片，为1时把孤立点合在一个分片其他分片不变
    std::set<std::set<std::string>> mergedSets = merge_high_similarity_sets(twoHopSets, merge_threshold);
    std::cout << "Number of merged sets: " << mergedSets.size() << std::endl;
    std::cout << "Number of twoHopSets: " << twoHopSets.size() << std::endl;
    
    int part = 0;
    // 通过合并后的种子集合生成分片
    std::ofstream tmpfile("./route/person2part.txt");
    for(auto i: mergedSets) {
        std::string outPath = "./output/";
        outPath+=std::to_string(part);
        std::cout<<"part"<<part<<std::endl;
        std::cout<<"twohop persons num:"<<i.size()<<std::endl;
        for(auto it=i.begin();it!=i.end();it++){
            // outPath += *it;
            // outPath += "_";
            tmpfile<<*it<<'|'<<part<<std::endl;
        }
        std::set<std::string> comment;
        std::set<std::string> post;
        getDivide(i, inputdir, outPath, comment, post, knows);
        std::cout<<"part persons num:"<<i.size()<<std::endl<<std::endl;
        part++;
    }

}

void partition3_1(int input_threshhold, int input_similarity, std::string inputdir, std::map<std::string, Person> people, std::map<std::string, std::map<std::string,std::string>> knows,std::map<std::string, std::set<std::string>> two_hop_neighbors, std::map<std::string, int> neighbor_counts) {

    // person2part还有问题！
    // 划分方法3.1 一次load：3.1贪心算法合并person，3.2合并相似度高的两跳person集
    std::unordered_map<std::string,std::unordered_set<std::string>> person_hasInterest_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> person_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_studyAt_organisation;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_workAt_organisation;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> forum_hasMember_person_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> forum_hasModerator_person_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> forum_hasModerator_person;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasCreator_person_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasCreator_person;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasCreator_person_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasCreator_person;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasTag_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasTag_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_set<std::string>> forum_containerOf_post_r;
    std::unordered_map<std::string, std::unordered_set<std::string>> reply1;
    std::unordered_map<std::string, std::unordered_set<std::string>> reply2;
    std::unordered_map<std::string, std::string> person_0_0;
    std::unordered_map<std::string, std::string> comment_0_0;
    std::unordered_map<std::string, std::string> post_0_0;
    std::unordered_map<std::string, std::string> forum_0_0;

    loadExknows(inputdir,knows, person_hasInterest_tag,person_isLocatedIn_place,person_studyAt_organisation,person_workAt_organisation,forum_hasMember_person_r,forum_hasModerator_person_r,forum_hasModerator_person,comment_hasCreator_person_r,comment_hasCreator_person,post_hasCreator_person_r,post_hasCreator_person,post_hasTag_tag,post_isLocatedIn_place,comment_hasTag_tag,comment_isLocatedIn_place,forum_containerOf_post_r,reply1,reply2,person_0_0,comment_0_0,post_0_0,forum_0_0);

    
    std::ofstream outfile("./route/message2person.txt");
    for(auto i:comment_hasCreator_person){
        for(auto j:comment_hasCreator_person[i.first]){
            outfile<<i.first<<'|'<<j<<std::endl;
        }
    }
    for(auto i:post_hasCreator_person){
        for(auto j:post_hasCreator_person[i.first]){
            outfile<<i.first<<'|'<<j<<std::endl;
        }
    }
    outfile.close();

    std::set<std::set<std::string>> twoHopSets;   // 记录所有分片的两步邻居

    std::set<std::string> peopleSet;    // 记录剩余未做种子分片的person，初始为所有person
    for(auto it=people.begin();it!=people.end();it++){
        peopleSet.insert(it->first);
    }

    int testCenter=0;
    int testtwohop=0;

    while(peopleSet.size()!=0) {    // 存在未作种子的person
        // 取seed，即一跳邻居最多的person
        std::string seed = getBiggestOneHop(peopleSet,neighbor_counts);
        if(seed=="") continue;
        
        //将seed及其两跳邻居加入分片
        std::set<std::string> twoHopResults;
        twoHopResults.insert(seed);
        twoHopResults.insert(two_hop_neighbors[seed].begin(), two_hop_neighbors[seed].end());
        

        // 将seed加入种子集即center
        std::set<std::string> center;
        center.insert(seed);
        // 获取seed的一跳邻居，即候选种子，不包含已经作为center的person
        std::set<std::string> cNeighbor = getKnowsNeighbor1(center,knows,peopleSet);
        for (const auto& person : center) {
            cNeighbor.erase(person);
        }

        int threshhold = twoHopResults.size();  // 阈值初始为seed的两跳邻居数目
        while(cNeighbor.size()!=0&&threshhold<input_threshhold) {

            // 在邻居中找到一跳邻居最小的person
            std::string smallestOneHop = getSmallestOneHop(cNeighbor,neighbor_counts);
            if (smallestOneHop == "") { // cNeighbor为空时返回为""，正常情况下循环直接结束其实不会执行
                break;
            }
            // 更新阈值
            threshhold = calculateSetTwoHopNeighborNumNew(threshhold, smallestOneHop, twoHopResults, two_hop_neighbors);

            // 符合阈值则加入center和两步邻居
            center.insert(smallestOneHop);
            twoHopResults.insert(smallestOneHop);
            twoHopResults.insert(two_hop_neighbors[smallestOneHop].begin(), two_hop_neighbors[smallestOneHop].end());
            // peopleSet.erase(smallestOneHop); //这里可删可不删，cNeighbor.erase(center);会删掉

            // 更新候选种子集，这里已经在两步邻居里的一步邻居还要做种子，相当于从比较小的扩张
            cNeighbor = getKnowsNeighbor1(center,knows,peopleSet);
            // 保证做过center的不会再被选为候选种子
            for (const auto& person : center) {
                cNeighbor.erase(person);
            }
            
        }
        
        // 插入到twoHopSets中待合并
        twoHopSets.insert(twoHopResults);
        
        // 保证做过种子的不会再在其他分片做种子
        for (const auto& person : center) {
            peopleSet.erase(person);
        }

        testCenter+=center.size();
        testtwohop+=twoHopResults.size();
    }
    std::cout<<"all center person num:"<<testCenter<<std::endl;
    std::cout<<"all twohop person num:"<<testtwohop<<std::endl;
    
    // 合并相似度高的两跳person集
    double merge_threshold = input_similarity;  // 相似度阈值,为0时合为一个分片，为1时把孤立点合在一个分片其他分片不变
    std::set<std::set<std::string>> mergedSets = merge_high_similarity_sets(twoHopSets, merge_threshold);
    std::cout << "Number of merged sets: " << mergedSets.size() << std::endl;
    std::cout << "Number of twoHopSets: " << twoHopSets.size() << std::endl;
    
    int part = 0;
    // 通过合并后的种子集合生成分片
    std::ofstream tmpfile("./route/person2part.txt");
    for(auto i: mergedSets) {
        std::string outPath = "./output/";
        outPath+=std::to_string(part);
        std::cout<<"part"<<part<<std::endl;
        std::cout<<"twohop persons num:"<<i.size()<<std::endl;
        for(auto it=i.begin();it!=i.end();it++){
            // outPath += *it;
            // outPath += "_";
            tmpfile<<*it<<'|'<<part<<std::endl;
        }
        std::set<std::string> comment;
        std::set<std::string> post;
        divideInMemory(i, outPath, comment, post, knows, person_hasInterest_tag, person_isLocatedIn_place, person_studyAt_organisation, person_workAt_organisation, forum_hasMember_person_r, forum_hasModerator_person_r, forum_hasModerator_person, comment_hasCreator_person_r, comment_hasCreator_person, post_hasCreator_person_r, post_hasCreator_person, post_hasTag_tag, post_isLocatedIn_place, comment_hasTag_tag, comment_isLocatedIn_place, forum_containerOf_post_r, reply1, reply2, person_0_0, comment_0_0, post_0_0, forum_0_0);
        std::cout<<"part persons num:"<<i.size()<<std::endl<<std::endl;
        part++;
    }
}

void partitionMultiThreads(int input_threshhold, int input_similarity, std::string inputdir, std::map<std::string, Person> people, std::map<std::string, std::map<std::string,std::string>> knows,std::map<std::string, std::set<std::string>> two_hop_neighbors, std::map<std::string, int> neighbor_counts) {

    // 划分方法3.1 一次load：3.1贪心算法合并person，3.2合并相似度高的两跳person集
    std::unordered_map<std::string,std::unordered_set<std::string>> person_hasInterest_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> person_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_studyAt_organisation;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_workAt_organisation;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> forum_hasMember_person_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> forum_hasModerator_person_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> forum_hasModerator_person;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasCreator_person_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasCreator_person;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasCreator_person_r;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasCreator_person;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasTag_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> post_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasTag_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> comment_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_set<std::string>> forum_containerOf_post_r;
    std::unordered_map<std::string, std::unordered_set<std::string>> reply1;
    std::unordered_map<std::string, std::unordered_set<std::string>> reply2;
    std::unordered_map<std::string, std::string> person_0_0;
    std::unordered_map<std::string, std::string> comment_0_0;
    std::unordered_map<std::string, std::string> post_0_0;
    std::unordered_map<std::string, std::string> forum_0_0;

    loadExknows(inputdir,knows, person_hasInterest_tag,person_isLocatedIn_place,person_studyAt_organisation,person_workAt_organisation,forum_hasMember_person_r,forum_hasModerator_person_r,forum_hasModerator_person,comment_hasCreator_person_r,comment_hasCreator_person,post_hasCreator_person_r,post_hasCreator_person,post_hasTag_tag,post_isLocatedIn_place,comment_hasTag_tag,comment_isLocatedIn_place,forum_containerOf_post_r,reply1,reply2,person_0_0,comment_0_0,post_0_0,forum_0_0);

    
    std::ofstream outfile("./route/message2person.txt");
    for(auto i:comment_hasCreator_person){
        for(auto j:comment_hasCreator_person[i.first]){
            outfile<<i.first<<'|'<<j<<std::endl;
        }
    }
    for(auto i:post_hasCreator_person){
        for(auto j:post_hasCreator_person[i.first]){
            outfile<<i.first<<'|'<<j<<std::endl;
        }
    }
    outfile.close();

    std::set<std::set<std::string>> twoHopSets;   // 记录所有分片的两步邻居

    std::set<std::string> peopleSet;    // 记录剩余未做种子分片的person，初始为所有person
    for(auto it=people.begin();it!=people.end();it++){
        peopleSet.insert(it->first);
    }

    int testCenter=0;
    int testtwohop=0;

    while(peopleSet.size()!=0) {    // 存在未作种子的person
        // 取seed，即一跳邻居最多的person
        std::string seed = getBiggestOneHop(peopleSet,neighbor_counts);
        if(seed=="") continue;
        
        //将seed及其两跳邻居加入分片
        std::set<std::string> twoHopResults;
        twoHopResults.insert(seed);
        twoHopResults.insert(two_hop_neighbors[seed].begin(), two_hop_neighbors[seed].end());
        

        // 将seed加入种子集即center
        std::set<std::string> center;
        center.insert(seed);
        // 获取seed的一跳邻居，即候选种子，不包含已经作为center的person
        std::set<std::string> cNeighbor = getKnowsNeighbor1(center,knows,peopleSet);
        for (const auto& person : center) {
            cNeighbor.erase(person);
        }

        int threshhold = twoHopResults.size();  // 阈值初始为seed的两跳邻居数目
        while(cNeighbor.size()!=0&&threshhold<input_threshhold) {

            // 在邻居中找到一跳邻居最小的person
            std::string smallestOneHop = getSmallestOneHop(cNeighbor,neighbor_counts);
            if (smallestOneHop == "") { // cNeighbor为空时返回为""，正常情况下循环直接结束其实不会执行
                break;
            }
            // 更新阈值
            threshhold = calculateSetTwoHopNeighborNumNew(threshhold, smallestOneHop, twoHopResults, two_hop_neighbors);

            // 符合阈值则加入center和两步邻居
            center.insert(smallestOneHop);
            twoHopResults.insert(smallestOneHop);
            twoHopResults.insert(two_hop_neighbors[smallestOneHop].begin(), two_hop_neighbors[smallestOneHop].end());
            // peopleSet.erase(smallestOneHop); //这里可删可不删，cNeighbor.erase(center);会删掉

            // 更新候选种子集，这里已经在两步邻居里的一步邻居还要做种子，相当于从比较小的扩张
            cNeighbor = getKnowsNeighbor1(center,knows,peopleSet);
            // 保证做过center的不会再被选为候选种子
            for (const auto& person : center) {
                cNeighbor.erase(person);
            }
            
        }
        
        // 插入到twoHopSets中待合并
        twoHopSets.insert(twoHopResults);
        
        // 保证做过种子的不会再在其他分片做种子
        for (const auto& person : center) {
            peopleSet.erase(person);
        }

        testCenter+=center.size();
        testtwohop+=twoHopResults.size();
    }
    std::cout<<"all center person num:"<<testCenter<<std::endl;
    std::cout<<"all twohop person num:"<<testtwohop<<std::endl;
    
    // 合并相似度高的两跳person集
    double merge_threshold = input_similarity;  // 相似度阈值,为0时合为一个分片，为1时把孤立点合在一个分片其他分片不变
    std::set<std::set<std::string>> mergedSets = merge_high_similarity_sets(twoHopSets, merge_threshold);
    std::cout << "Number of merged sets: " << mergedSets.size() << std::endl;
    std::cout << "Number of twoHopSets: " << twoHopSets.size() << std::endl;
    
    int part = 0;
    // 通过合并后的种子集合生成分片
    std::ofstream tmpfile("./route/person2part.txt");
    for(auto i: mergedSets) {
        std::string outPath = "./output/";
        outPath+=std::to_string(part);
        std::cout<<"part"<<part<<std::endl;
        std::cout<<"twohop persons num:"<<i.size()<<std::endl;
        for(auto it=i.begin();it!=i.end();it++){
            // outPath += *it;
            // outPath += "_";
            tmpfile<<*it<<'|'<<part<<std::endl;
        }
        std::set<std::string> comment;
        std::set<std::string> post;
        divideInMemory(i, outPath, comment, post, knows, person_hasInterest_tag, person_isLocatedIn_place, person_studyAt_organisation, person_workAt_organisation, forum_hasMember_person_r, forum_hasModerator_person_r, forum_hasModerator_person, comment_hasCreator_person_r, comment_hasCreator_person, post_hasCreator_person_r, post_hasCreator_person, post_hasTag_tag, post_isLocatedIn_place, comment_hasTag_tag, comment_isLocatedIn_place, forum_containerOf_post_r, reply1, reply2, person_0_0, comment_0_0, post_0_0, forum_0_0);
        std::cout<<"part persons num:"<<i.size()<<std::endl<<std::endl;
        part++;
    }
}

// 运行长度编码压缩函数
std::string compressString(const std::string& str) {
    std::string compressed;
    int count = 1;

    for (size_t i = 1; i <= str.length(); ++i) {
        if (i < str.length() && str[i] == str[i - 1]) {
            ++count;
        } else {
            compressed += str[i - 1];
            compressed += std::to_string(count);
            count = 1;
        }
    }

    return compressed.length() < str.length() ? compressed : str;
}

// 运行长度编码解压缩函数
std::string decompressString(const std::string& str) {
    std::string decompressed;
    for (size_t i = 0; i < str.length(); ++i) {
        char ch = str[i];
        std::string countStr;
        while (i + 1 < str.length() && isdigit(str[i + 1])) {
            countStr += str[++i];
        }
        int count = std::stoi(countStr);
        decompressed.append(count, ch);
    }
    return decompressed;
}

void loadAll(std::string filename,
                std::unordered_map<std::string,std::unordered_set<std::string>>& person_hasInterest_tag,
                std::unordered_map<std::string,std::unordered_set<std::string>>& person_isLocatedIn_place,
                std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& person_studyAt_organisation,
                std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& person_workAt_organisation,
                std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& person_likes_comment_r,
                std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& person_likes_post_r,
                // std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& forum_hasMember_person_r,
                // std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& forum_hasMember_person, //new
                // std::unordered_map<std::string,std::unordered_set<std::string>>& forum_hasModerator_person_r,
                // std::unordered_map<std::string,std::string>& forum_hasModerator_person,
                std::unordered_map<std::string,std::unordered_set<std::string>>& comment_hasCreator_person_r,
                std::unordered_map<std::string,std::string>& comment_hasCreator_person,
                std::unordered_map<std::string,std::unordered_set<std::string>>& post_hasCreator_person_r,
                std::unordered_map<std::string,std::string>& post_hasCreator_person,
                std::unordered_map<std::string,std::unordered_set<std::string>>& post_hasTag_tag,
                std::unordered_map<std::string,std::unordered_set<std::string>>& post_isLocatedIn_place,
                std::unordered_map<std::string,std::unordered_set<std::string>>& comment_hasTag_tag,
                std::unordered_map<std::string,std::unordered_set<std::string>>& comment_isLocatedIn_place,
                std::unordered_map<std::string,std::unordered_set<std::string>>& forum_containerOf_post, //new
                std::unordered_map<std::string,std::unordered_set<std::string>>& forum_containerOf_post_r, //new
                std::unordered_map<std::string,std::unordered_set<std::string>>& comment_replyOf_comment_r,
                std::unordered_map<std::string,std::unordered_set<std::string>>& comment_replyOf_post_r,
                std::unordered_map<std::string,std::unordered_set<std::string>>& comment_replyOf_comment,
                std::unordered_map<std::string,std::unordered_set<std::string>>& comment_replyOf_post,
                std::unordered_map<std::string, std::string>& person_0_0,
                std::unordered_map<std::string, std::string>& comment_0_0,
                std::unordered_map<std::string, std::string>& post_0_0,
                std::unordered_map<std::string, std::string>& forum_0_0) {

    // 获取文件夹路径
    // std::string directory = filename.substr(0, filename.find_last_of('/'));
    std::string directory = filename;

    std::cout<<"inputdir:"<<directory<<std::endl;

    // std::ifstream file(directory+"person_0_0.csv");
    // std::string line;

    // while (std::getline(file, line)) {
    //     std::stringstream ss(line);
    //     std::string item;
    //     std::vector<std::string> items;

    //     while (std::getline(ss, item, '|')) {
    //         items.push_back(item);
    //     }

    //     if (items.size() != 10) {
    //         std::cerr << "Invalid line: " << line << std::endl;
    //         continue;
    //     }
    //     person[items[0]]=line;
    // }
    
    // person_hasInterest_tag
    std::ifstream file(directory+"/person_hasInterest_tag_0_0.csv");
    std::string line;

    // std::getline(file, line); // Skip header line
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
        person_hasInterest_tag[items[0]].insert(items[1]);
    }
    file.close();
    
    // person_isLocatedIn_city
    
    file.open(directory+"/person_isLocatedIn_place_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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
        person_isLocatedIn_place[items[0]].insert(items[1]);
    }
    file.close();
    
    // person_studyAt_university
    
    file.open(directory+"/person_studyAt_organisation_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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
        person_studyAt_organisation[items[0]][items[1]] = items[2];
    }
    file.close();
    
    
    // person_workAt_company
    
    file.open(directory+"/person_workAt_organisation_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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
        person_workAt_organisation[items[0]][items[1]] = items[2];
    }
    file.close();

    // person_likes_comment
    
    file.open(directory+"/person_likes_comment_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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
    
    // std::getline(file, line); // Skip header line
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
    
    
    // // forum_hasMember_person
    
    // file.open(directory+"/forum_hasMember_person_0_0.csv");
    
    // // std::getline(file, line); // Skip header line
    // while (std::getline(file, line)) {
    //     std::stringstream ss(line);
    //     std::string item;
    //     std::vector<std::string> items;

    //     while (std::getline(ss, item, '|')) {
    //         items.push_back(item);
    //     }

    //     if (items.size() != 3) {
    //         std::cerr << "Invalid line: " << line << std::endl;
    //         continue;
    //     }
    //     forum_hasMember_person_r[items[1]][items[0]] = items[2];
    //     forum_hasMember_person[items[0]][items[1]] = items[2];
    // }
    // file.close();
    
    // // forum_hasModerator_person

    // file.open(directory+"/forum_hasModerator_person_0_0.csv");
    
    // // std::getline(file, line); // Skip header line
    // while (std::getline(file, line)) {
    //     std::stringstream ss(line);
    //     std::string item;
    //     std::vector<std::string> items;

    //     while (std::getline(ss, item, '|')) {
    //         items.push_back(item);
    //     }

    //     if (items.size() != 2) {
    //         std::cerr << "Invalid line: " << line << std::endl;
    //         continue;
    //     }
    //     forum_hasModerator_person_r[items[1]].insert(items[0]);
    //     forum_hasModerator_person[items[0]]=items[1];
    // }
    // file.close();

    // comment hasCreator person
    file.open(directory+"/comment_hasCreator_person_0_0.csv");
    // file2<<"Comment.id|Person.id"<<std::endl;    //为什么没有输出到呢，open会覆盖

    // std::getline(file, line); // Skip header line
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
    // file2<<"Post.id|Person.id"<<std::endl;

    // std::getline(file, line); // Skip header line
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
    
    
	// post_hasTag_tag
    file.open(directory+"/post_hasTag_tag_0_0.csv");
    // std::set <std::string> newPerson;
    
    // std::getline(file, line); // Skip header line
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

	// post_isLocatedIn_place
    file.open(directory+"/post_isLocatedIn_place_0_0.csv");
    // std::set <std::string> newPerson;
    
    // std::getline(file, line); // Skip header line
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
    
    // std::getline(file, line); // Skip header line
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
    
    // std::getline(file, line); // Skip header line
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
    
    
	// forum_containerOf_post
    file.open(directory+"/forum_containerOf_post_0_0.csv");
    std::set <std::string> newForum;
    
    // std::getline(file, line); // Skip header line
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
        forum_containerOf_post[items[0]].insert(items[1]);
        forum_containerOf_post_r[items[1]].insert(items[0]);
    }
    file.close();

	// forum
    file.open(directory+"/forum_0_0.csv");
    
    // std::getline(file, line); // Skip header line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }
        forum_0_0[items[0]] = line;
    }
    file.close();

	// person
    file.open(directory+"/person_0_0.csv");
    
    // std::getline(file, line); // Skip header line
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

	// comment
    file.open(directory+"/comment_0_0.csv");
    
    // std::getline(file, line); // Skip header line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }
        comment_0_0[items[0]] = compressString(line);
    }
    file.close();

	// post
    file.open(directory+"/post_0_0.csv");
    
    // std::getline(file, line); // Skip header line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }
        post_0_0[items[0]] = compressString(line);
    }
    file.close();
    std::cout<<"load finished"<<std::endl;

}

void divideMessage(int part, std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasCreator_person_r, 
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasCreator_person_r,
                                                    std::unordered_map<std::string,std::string> comment_hasCreator_person, 
                                                    std::unordered_map<std::string,std::string> post_hasCreator_person,
                                                    std::unordered_map<std::string,std::string> comment_0_0,
                                                    std::unordered_map<std::string,std::string> post_0_0,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> comment_replyOf_comment_r,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> comment_replyOf_post_r,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> comment_replyOf_comment,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> comment_replyOf_post,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasTag_tag,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasTag_tag,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> comment_isLocatedIn_place,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> post_isLocatedIn_place) {
    // 将message按照personid%part划分到part集合中
    std::vector<std::unordered_set<std::string>> commentSetArray(part);
    std::vector<std::unordered_set<std::string>> postSetArray(part);
    // 对comment和post按照personid%part划分
    for(auto i:comment_hasCreator_person_r){
        long long int index = std::stol(i.first)%part;
        commentSetArray[index].insert(i.second.begin(),i.second.end());
    }
    for(auto i:post_hasCreator_person_r){
        long long int index = std::stol(i.first)%part;
        postSetArray[index].insert(i.second.begin(),i.second.end());
    }
    
    std::string outPath="./output/";
    if (!createDirectory("./output/")) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    if (!createDirectory("./route/")) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    std::ofstream message2partf("./route/message2part.csv");
    std::unordered_map<std::string,std::string> message2part;
    for(int i=0;i<part;i++){
        for(auto j:commentSetArray[i]){
            message2part[j]=std::to_string(i);
        }
        for(auto j:postSetArray[i]){
            message2part[j]=std::to_string(i);
        }
    }
    for(auto i:message2part){
        message2partf<<i.first<<'|'<<i.second<<std::endl;
    }
    message2partf.close();
    for(int i=0;i<part;i++){
        outPath="./output/";
        outPath+=std::to_string(i);
        // 如果文件夹不存在，则创建它
        if (!createDirectory(outPath)) {
            // std::cerr << "Failed or no need to create directory." << std::endl;
        }
        std::cout<<outPath<<std::endl;
        std::ofstream outfile(outPath+"/comment_hasCreator_person_0_0.csv");
        std::ofstream outfile0(outPath+"/post_hasCreator_person_0_0.csv");
        std::ofstream outfile1(outPath+"/comment_0_0.csv");
        std::ofstream outfile2(outPath+"/post_0_0.csv");

        std::ofstream file3_1(outPath+"/comment_replyOf_comment_0_0.csv");
        std::ofstream file3_2(outPath+"/comment_replyOf_post_0_0.csv");

        
        std::ofstream outfile4(outPath+"/comment_isLocatedIn_place_0_0.csv");
        std::ofstream outfile5(outPath+"/post_isLocatedIn_place_0_0.csv");
        std::ofstream outfile6(outPath+"/comment_hasTag_tag_0_0.csv");
        std::ofstream outfile7(outPath+"/post_hasTag_tag_0_0.csv");

        if (outfile.is_open()) {
            std::cout << "文件成功打开，现在可以写入数据。" << std::endl;
            // 在这里写入数据到文件
        } else {
            std::cout << "文件打开失败，请检查路径和权限：" << outPath+"/post_hasCreator_person_0_0.csv" << std::endl;
        }
        // 去重comment2comment
        std::unordered_map<std::string,std::unordered_map<std::string,int>> visited; 
        // 去重comment2post
        std::unordered_map<std::string,std::unordered_map<std::string,int>> visited1; 
        
        std::unordered_set<std::string> tmpComment;
        std::unordered_set<std::string> tmpPost;
        std::unordered_set<std::string> tmpPerson;
        tmpComment.insert(commentSetArray[i].begin(),commentSetArray[i].end());
        tmpPost.insert(postSetArray[i].begin(),postSetArray[i].end());
        for(auto k:commentSetArray[i]) {
            // comment的前向comment和边
            tmpComment.insert(comment_replyOf_comment_r[k].begin(),comment_replyOf_comment_r[k].end());
            for(auto kk:comment_replyOf_comment_r[k]){
                if(visited[kk].find(k)==visited[kk].end()) {
                        visited[kk][k] = 1;
                        file3_1<<kk<<'|'<<k<<std::endl;
                }
            }
            // comment的后向comment和边
            tmpComment.insert(comment_replyOf_comment[k].begin(),comment_replyOf_comment[k].end());
            for(auto kk:comment_replyOf_comment[k]){
                if(visited[k].find(kk)==visited[k].end()) {
                        visited[k][kk] = 1;
                        file3_1<<k<<'|'<<kk<<std::endl;
                }
            }
            // comment的后向post和边
            tmpPost.insert(comment_replyOf_post[k].begin(),comment_replyOf_post[k].end());
            for(auto kk:comment_replyOf_post[k]){
                if(visited1[k].find(kk)==visited1[k].end()) {
                        visited1[k][kk] = 1;
                        file3_2<<k<<'|'<<kk<<std::endl;
                }
            }
        }
        for(auto k:postSetArray[i]) {
            // post的前向
            tmpComment.insert(comment_replyOf_post_r[k].begin(),comment_replyOf_post_r[k].end());
            for(auto kk:comment_replyOf_post_r[k]){
                if(visited1[kk].find(k)==visited1[kk].end()) {
                        visited1[kk][k] = 1;
                        file3_2<<kk<<'|'<<k<<std::endl;
                }
            }
        }
        
        std::unordered_map<std::string,std::vector<std::string>> commentPath;
        for (auto k:commentSetArray[i]) {
            std::vector<std::string> tmp;
            dfs(k, comment_replyOf_comment, tmp, commentPath);
        }
        for(auto k=commentPath.begin();k!=commentPath.end();k++){
            std::string from = k->first;
            std::vector<std::string> tmp = k->second;
            if(tmp.size()==1) {
                // comment没有comment邻居
                // ！！！这里也要存post
                tmpComment.insert(tmp[0]); // 里面本来就存了
                // file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                for(std::string j : comment_replyOf_comment[tmp[0]]) {
                    tmpPost.insert(j);
                    if(visited1[tmp[0]].find(j)==visited1[tmp[0]].end()) {
                        visited1[tmp[0]][j] = 1;
                        file3_2<<tmp[0]<<'|'<<j<<std::endl;
                    }
                    // file3_2<<tmp[0]<<'|'<<j<<std::endl;
                }
            } else if(tmp.size()==2) {
                // comment只有一步comment
                tmpComment.insert(tmp[1]);
                if(visited[tmp[0]].find(tmp[1])==visited[tmp[0]].end()) {
                    visited[tmp[0]][tmp[1]] = 1;
                    file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                }
                // file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                for(std::string j : comment_replyOf_post[tmp[1]]) {
                    // comment存在comment链
                    tmpPost.insert(j);
                    if(visited1[tmp[1]].find(j)==visited1[tmp[1]].end()) {
                        visited1[tmp[1]][j] = 1;
                        file3_2<<tmp[1]<<'|'<<j<<std::endl;
                    }
                    // file3_2<<tmp[1]<<'|'<<j<<std::endl;
                }
            } else {
                // 顶点有多步邻居
                tmpComment.insert(tmp[1]);
                if(visited[tmp[0]].find(tmp[1])==visited[tmp[0]].end()) {
                    visited[tmp[0]][tmp[1]] = 1;
                    file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                }
                // 这里把末端邻居合并为post
                for(std::string j : comment_replyOf_post[tmp[2]]) {
                    tmpPost.insert(j);
                    if(visited1[tmp[1]].find(j)==visited1[tmp[1]].end()) {
                        visited1[tmp[1]][j] = 1;
                        file3_2<<tmp[1]<<'|'<<j<<std::endl;
                    }
                    // 压缩了路径，不过好像导致comment2post变多了？
                }
            }
        }

        // 写入comment
        for(auto k:tmpComment){
            outfile<<k<<'|'<<comment_hasCreator_person[k]<<std::endl;
            outfile1<<comment_0_0[k]<<std::endl;
            for(auto j:comment_isLocatedIn_place[k]){
                outfile4<<k<<'|'<<j<<std::endl;
            }
            for(auto j:comment_hasTag_tag[k]){
                outfile6<<k<<'|'<<j<<std::endl;
            }
        }
        // 写入post
        for(auto k:tmpPost){
            outfile0<<k<<'|'<<post_hasCreator_person[k]<<std::endl;
            outfile2<<post_0_0[k]<<std::endl;
            for(auto j:post_isLocatedIn_place[k]){
                outfile5<<k<<'|'<<j<<std::endl;
            }
            for(auto j:post_hasTag_tag[k]){
                outfile7<<k<<'|'<<j<<std::endl;
            }
        }
        
        outfile.close();
        outfile0.close();
        outfile1.close();
        outfile2.close();
        file3_1.close();
        file3_2.close();
        outfile4.close();
        outfile5.close();
        outfile6.close();
        outfile7.close();

    }

}

// 合并两个 map，用于合并comment和post的replyOf关系
std::unordered_map<std::string, std::unordered_set<std::string>> merge_maps(
    const std::unordered_map<std::string, std::unordered_set<std::string>> &map1,
    const std::unordered_map<std::string, std::unordered_set<std::string>> &map2) {
    
    std::unordered_map<std::string, std::unordered_set<std::string>> result;

    // 将 map1 的所有键值对插入到 result 中
    for (const auto &pair : map1) {
        result[pair.first] = pair.second;
    }

    // 遍历 map2，合并到 result 中
    for (const auto &pair : map2) {
        const std::string &key = pair.first;
        const std::unordered_set<std::string> &set2 = pair.second;
        
        if (result.find(key) != result.end()) {
            // 如果 result 中已经有这个 key，则合并 unordered_set
            std::unordered_set<std::string> &set1 = result[key];
            set1.insert(set2.begin(), set2.end());
        } else {
            // 如果 result 中没有这个 key，则直接插入
            result[key] = set2;
        }
    }

    return result;
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
     // 获取文件夹路径
    // std::string directory = filename.substr(0, filename.find_last_of('/'));
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
    
    // std::getline(file, line); // Skip header line
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
    
    // std::getline(file, line); // Skip header line
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

void loadMessageCreationDate(std::string filename, std::unordered_map<std::string, std::string>& comment_0_0, std::unordered_map<std::string, std::string>& post_0_0) {
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
        comment_0_0[items[0]] = items[1];
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
        post_0_0[items[0]] = items[2];
    }
    file.close();
}

void loadMessageProperty(std::string filename, 
                 std::unordered_map<std::string, std::string>& comment_0_0, 
                 std::unordered_map<std::string, std::string>& post_0_0,
                 std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& person_likes_comment_r,
                 std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& person_likes_post_r,
                //  std::unordered_map<std::string,std::unordered_set<std::string>>& post_hasTag_tag,
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

    file.open(directory+"/person_likes_comment_0_0.csv");
    // std::getline(file, line); // Skip header line
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
    
    // std::getline(file, line); // Skip header line
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

    // // post_hasTag_tag
    // file.open(directory+"/post_hasTag_tag_0_0.csv");
    // // std::set <std::string> newPerson;
    
    // // std::getline(file, line); // Skip header line
    // while (std::getline(file, line)) {
    //     std::stringstream ss(line);
    //     std::string item;
    //     std::vector<std::string> items;

    //     while (std::getline(ss, item, '|')) {
    //         items.push_back(item);
    //     }

    //     if (items.size() != 2) {
    //         std::cerr << "Invalid line: " << line << std::endl;
    //         continue;
    //     }
        
    //     post_hasTag_tag[items[0]].insert(items[1]);
    // }
    // file.close();

	// post_isLocatedIn_place
    file.open(directory+"/post_isLocatedIn_place_0_0.csv");
    // std::set <std::string> newPerson;
    
    // std::getline(file, line); // Skip header line
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
    
    // std::getline(file, line); // Skip header line
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
    
    // std::getline(file, line); // Skip header line
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
    std::cout<<"load finished"<<std::endl;
  
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
    // 便于计算先合并comment和post的replyOf
    // std::unordered_map<std::string,std::unordered_set<std::string>> replySet = merge_maps(comment_replyOf_comment,comment_replyOf_post);
    // 处理每两个person点对
    for(auto i:person_0_0) {
        // 初始化person i的comment集
        std::unordered_set<std::string> iSet = comment_hasCreator_person_r[i.first];
        //其实这里不用,post不会被回复
        // iSet.insert(post_hasCreator_person_r[i.first].begin(),post_hasCreator_person_r[i.first].end());  
        // dfs得到所有replyOf*邻居,一个message只会reply一个，所以直接统计所有replyOf*邻居
        std::set<std::string> i_reply1,i_reply2;
        for(auto ii: iSet) {
            std::unordered_map<std::string,std::unordered_map<std::string,int>> l_visited,l_visited1; 
            // dfs_c(ii,replySet,i_reply,l_visited);   //拆分成两个应该也一样
            dfs_c(ii,comment_replyOf_comment,i_reply1,l_visited);
            dfs_c(ii,comment_replyOf_post,i_reply2,l_visited1);
        }
        for(auto j: i_reply1) {
            // comment
            std::string to = comment_hasCreator_person[j];
            if(i.first!=to){
                g.addEdge(i.first,to,0.5);
            }
        }
        for(auto j: i_reply2) {
            // post
            std::string to = post_hasCreator_person[j];
            if(i.first!=to){
                g.addEdge(i.first,to,0.5);
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    
    std::cout<<"create louvain graph finished,uses "<<duration.count()<<"ms"<<std::endl;
    
}

void divideCommunity(int part,int messageNum,Graph& g,
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
    // GetCurrentProcessMemoryUsage();
    auto start = std::chrono::high_resolution_clock::now();
    std::unordered_set<std::string> communityMessageSet; //存储当前社区中person发布的message，用于判断阈值
    // unordered_map<string, std::unordered_set<std::string>> communitiesSet;  //存储每个社区中person
    std::unordered_set<std::string> pPerson;  //存储已分配到社区中person
    int communityNo=0;
    int thresh = 0;
    // 分配社区中的person
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
        // communityMessageSet.insert(person_message_set[center].begin(),person_message_set[center].end());
        communityMessageSet.insert(comment_hasCreator_person_r[center].begin(),comment_hasCreator_person_r[center].end());
        communityMessageSet.insert(post_hasCreator_person_r[center].begin(),post_hasCreator_person_r[center].end());
        thresh = communityMessageSet.size();
        if(thresh>=(messageNum/part)) {
            std::cout<<"divide community phase1 finished,community size:"<<communitiesSet[to_string(communityNo)].size()<<std::endl;
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
            // communityMessageSet.insert(person_message_set[tmp].begin(),person_message_set[tmp].end());
            communityMessageSet.insert(comment_hasCreator_person_r[tmp].begin(),comment_hasCreator_person_r[tmp].end());
            communityMessageSet.insert(post_hasCreator_person_r[tmp].begin(),post_hasCreator_person_r[tmp].end());
            thresh = communityMessageSet.size();
            if(thresh>=(messageNum/part)) {
                std::cout<<"divide community phase1 finished,community size:"<<communitiesSet[to_string(communityNo)].size()<<std::endl;
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
                std::cout<<"divide community phase2 finished,community size:"<<communitiesSet[to_string(communityNo)].size()<<std::endl;
                communityNo++;
                communityMessageSet.clear();
                thresh = 0;
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout<<"divide community finished,uses "<<duration.count()<<"s"<<endl;
    // 将message按照personid%part划分到part集合中
    // std::vector<std::unordered_set<std::string>> commentSetArray(part);
    // std::vector<std::unordered_set<std::string>> postSetArray(part);
    std::string outPath="./output";
    std::string outPath1="./route";
    if (!createDirectory(outPath)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    outPath="./output/messages";
    if (!createDirectory(outPath)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    if (!createDirectory(outPath1)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    std::ofstream message2partf("./route/message2part.csv");
    std::ofstream person2partf("./route/person2messagePart.csv");
    int commuCnt = 0;
    // int mc_num = part/communityNo;// 每个社区的message分片数量
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
            // GetCurrentProcessMemoryUsage();
            outPath="./output/messages/";
            outPath+=std::to_string(i);
            // 如果文件夹不存在，则创建它
            if (!createDirectory(outPath)) {
                std::cerr << "Failed or no need to create directory."<<outPath << std::endl;
            }
            std::cout<<outPath<<std::endl;
            std::ofstream outfile(outPath+"/comment_hasCreator_person_0_0.csv");
            std::ofstream outfile0(outPath+"/post_hasCreator_person_0_0.csv");
            // std::ofstream outfile1(outPath+"/comment_0_0.csv");
            std::ofstream outfile1_1(outPath+"/comment_0_0_added.csv");
            // std::ofstream outfile2(outPath+"/post_0_0.csv");
            std::ofstream outfile2_1(outPath+"/post_0_0_added.csv");

            std::ofstream file3_1(outPath+"/comment_replyOf_comment_0_0.csv");
            std::ofstream file3_2(outPath+"/comment_replyOf_post_0_0.csv");

            
            // std::ofstream outfile4(outPath+"/comment_isLocatedIn_place_0_0.csv");
            // std::ofstream outfile5(outPath+"/post_isLocatedIn_place_0_0.csv");
            // std::ofstream outfile6(outPath+"/comment_hasTag_tag_0_0.csv");
            std::ofstream outfile7(outPath+"/post_hasTag_tag_0_0.csv");

            std::ofstream outfile8(outPath+"/forum_containerOf_post_0_0.csv");
            std::ofstream outfile9(outPath+"/forum_0_0.csv");

            std::ofstream outfile10(outPath+"/person_0_0.csv");

            // std::ofstream likes1(outPath+"/person_likes_comment_0_0.csv");
            // std::ofstream likes2(outPath+"/person_likes_post_0_0.csv");

            // std::ofstream person_comment(outPath+"/person-comment_0_0.csv");
            // std::ofstream person_post(outPath+"/person-post_0_0.csv");
            // std::cout<<"open"<<std::endl;
            // GetCurrentProcessMemoryUsage();

            if (outfile.is_open()) {
                std::cout << "文件成功打开，现在可以写入数据。" << std::endl;
                // 在这里写入数据到文件
            } else {
                std::cout << "文件打开失败，请检查路径和权限：" << outPath+"/post_hasCreator_person_0_0.csv" << std::endl;
            }
            // 写入person
            for(auto i:communitiesSet[std::to_string(i)]) {
                outfile10<<i<<std::endl;
            }
            outfile10.close();

            
            // 去重comment2comment
            std::unordered_map<std::string,std::unordered_map<std::string,int>> visited; 
            // 去重comment2post
            std::unordered_map<std::string,std::unordered_map<std::string,int>> visited1; 
            // 去重当前分片的comment、post、forum
            std::unordered_set<std::string> tmpComment;
            std::unordered_set<std::string> tmpPost;
            std::unordered_set<std::string> tmpForum;
            // 存储commentPath
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
                            file3_1<<kk<<'|'<<k<<std::endl;
                    }
                }
                // comment的后向comment和边
                tmpComment.insert(comment_replyOf_comment[k].begin(),comment_replyOf_comment[k].end());
                for(auto kk:comment_replyOf_comment[k]){
                    if(visited[k].find(kk)==visited[k].end()) {
                            visited[k][kk] = 1;
                            file3_1<<k<<'|'<<kk<<std::endl;
                    }
                }
                // comment的后向post和边
                tmpPost.insert(comment_replyOf_post[k].begin(),comment_replyOf_post[k].end());
                for(auto kk:comment_replyOf_post[k]){
                    if(visited1[k].find(kk)==visited1[k].end()) {
                            visited1[k][kk] = 1;
                            file3_2<<k<<'|'<<kk<<std::endl;
                    }
                }
            }
            for(auto k:postSetArray[i]) {
                // post的前向
                tmpComment.insert(comment_replyOf_post_r[k].begin(),comment_replyOf_post_r[k].end());
                for(auto kk:comment_replyOf_post_r[k]){
                    if(visited1[kk].find(k)==visited1[kk].end()) {
                            visited1[kk][k] = 1;
                            file3_2<<kk<<'|'<<k<<std::endl;
                    }
                }
            }
            
            // std::unordered_map<std::string,std::vector<std::string>> commentPath;
            for (auto k:commentSetArray[i]) {
                std::vector<std::string> tmp;
                dfs(k, comment_replyOf_comment, tmp, commentPath);
            }
            for(auto k=commentPath.begin();k!=commentPath.end();k++){
                std::string from = k->first;
                std::vector<std::string> tmp = k->second;
                if(tmp.size()==1) {
                    // comment没有comment邻居
                    // ！！！这里也要存post
                    tmpComment.insert(tmp[0]); // 里面本来就存了
                    // file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                    for(std::string j : comment_replyOf_post[tmp[0]]) {
                        tmpPost.insert(j);
                        if(visited1[tmp[0]].find(j)==visited1[tmp[0]].end()) {
                            visited1[tmp[0]][j] = 1;
                            file3_2<<tmp[0]<<'|'<<j<<std::endl;
                        }
                        // file3_2<<tmp[0]<<'|'<<j<<std::endl;
                    }
                } else if(tmp.size()==2) {
                    // comment只有一步comment
                    tmpComment.insert(tmp[1]);
                    if(visited[tmp[0]].find(tmp[1])==visited[tmp[0]].end()) {
                        visited[tmp[0]][tmp[1]] = 1;
                        file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                    }
                    // file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                    for(std::string j : comment_replyOf_post[tmp[1]]) {
                        // comment存在comment链
                        tmpPost.insert(j);
                        if(visited1[tmp[1]].find(j)==visited1[tmp[1]].end()) {
                            visited1[tmp[1]][j] = 1;
                            file3_2<<tmp[1]<<'|'<<j<<std::endl;
                        }
                        // file3_2<<tmp[1]<<'|'<<j<<std::endl;
                    }
                } else {
                    // 顶点有多步邻居
                    tmpComment.insert(tmp[1]);
                    if(visited[tmp[0]].find(tmp[1])==visited[tmp[0]].end()) {
                        visited[tmp[0]][tmp[1]] = 1;
                        file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                    }
                    // 这里把末端邻居合并为post
                    for(std::string j : comment_replyOf_post[tmp[2]]) {
                        tmpPost.insert(j);
                        if(visited1[tmp[1]].find(j)==visited1[tmp[1]].end()) {
                            visited1[tmp[1]][j] = 1;
                            file3_2<<tmp[1]<<'|'<<j<<std::endl;
                        }
                        // 压缩了路径，不过好像导致comment2post变多了？
                    }
                }
            }
            
            // std::cout<<"message extend finished"<<std::endl;
            // GetCurrentProcessMemoryUsage();

            // 写入comment
            for(auto k:tmpComment){
                // 这里改成中心comment存整行，非中心comment存id（需要加空列）
                // outfile1原版本
                // outfile1<<comment_0_0[k]<<std::endl;
                // outfile1新版本
                if(commentSetArray[i].find(k)!=commentSetArray[i].end()) {
                    outfile<<k<<'|'<<comment_hasCreator_person[k]<<std::endl;
                    // 是中心message
                    // person_comment<<k<<std::endl;
                } else {
                    outfile1_1<<k<<"|"<<comment_creationDate[k]<<std::endl;
                }
            }
            // 写入post
            for(auto k:tmpPost){
                // outfile2<<post_0_0[k]<<std::endl;
                if(postSetArray[i].find(k)!=postSetArray[i].end()) {
                    outfile0<<k<<'|'<<post_hasCreator_person[k]<<std::endl;
                    // 是中心message
                    // person_post<<k<<std::endl;
                } else {
                    outfile2_1<<k<<"|"<<post_creationDate[k]<<std::endl;
                }
                for(auto j:post_hasTag_tag[k]){
                    outfile7<<k<<'|'<<j<<std::endl;
                }
                for(auto j:post_containerOf_forum[k]){
                    outfile8<<j<<'|'<<k<<std::endl;
                    tmpForum.insert(j);
                }
            }
            
            // 写入forum
            for(auto k:tmpForum){
                outfile9<<k<<std::endl;
            }
            
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
            
            outfile.close();
            outfile0.close();
            // outfile1.close();
            outfile1_1.close();
            // outfile2.close();
            outfile2_1.close();
            file3_1.close();
            file3_2.close();
            // outfile4.close();
            // outfile5.close();
            // outfile6.close();
            outfile7.close();
            outfile8.close();
            outfile9.close();
            // likes1.close();
            // likes2.close();
            // person_comment.close();
            // person_post.close();
            // 强制释放内存
            // std::cout<<"close"<<std::endl;
            // GetCurrentProcessMemoryUsage();
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
        // GetCurrentProcessMemoryUsage();
        outPath="./output/messages/";
        outPath+=std::to_string(i);
        // // 如果文件夹不存在，则创建它
        // if (!createDirectory(outPath)) {
        //     std::cerr << "Failed or no need to create directory."<<outPath << std::endl;
        // }
        std::cout<<outPath<<std::endl;
        std::ofstream outfile1(outPath+"/comment_0_0.csv");
        std::ofstream outfile2(outPath+"/post_0_0.csv");
        
        std::ofstream outfile4(outPath+"/comment_isLocatedIn_place_0_0.csv");
        std::ofstream outfile5(outPath+"/post_isLocatedIn_place_0_0.csv");
        std::ofstream outfile6(outPath+"/comment_hasTag_tag_0_0.csv");
        // std::ofstream outfile7(outPath+"/post_hasTag_tag_0_0.csv");

        std::ofstream likes1(outPath+"/person_likes_comment_0_0.csv");
        std::ofstream likes2(outPath+"/person_likes_post_0_0.csv");

        std::ofstream person_comment(outPath+"/person-comment_0_0.csv");
        std::ofstream person_post(outPath+"/person-post_0_0.csv");
        
        

        // 开始划分message
        for(auto k:commentSetArray[i]) {
            outfile1<<comment_0_0[k]<<std::endl;
            for(auto j:comment_isLocatedIn_place[k]){
                outfile4<<k<<'|'<<j<<std::endl;
            }
            for(auto j:comment_hasTag_tag[k]){
                outfile6<<k<<'|'<<j<<std::endl;
            }
            for(auto j:person_likes_comment_r[k]){
                likes1<<j.first<<'|'<<k<<'|'<<j.second<<std::endl;
            }
            person_comment<<k<<std::endl;
        }
        for(auto k:postSetArray[i]) {
            outfile2<<post_0_0[k]<<std::endl;
            for(auto j:post_isLocatedIn_place[k]){
                outfile5<<k<<'|'<<j<<std::endl;
            }
            // for(auto j:post_hasTag_tag[k]){
            //     outfile7<<k<<'|'<<j<<std::endl;
            // }
            for(auto j:person_likes_post_r[k]){
                likes2<<j.first<<'|'<<k<<'|'<<j.second<<std::endl;
            }
            person_post<<k<<std::endl;
        }
        
        outfile1.close();
        outfile2.close();
        outfile4.close();
        outfile5.close();
        outfile6.close();
        // outfile7.close();

        likes1.close();
        likes2.close();
        person_comment.close();
        person_post.close();
        
    }
}

void divideMessageCommunity(int part, int messageNum, std::unordered_map<std::string,std::string> &person_0_0, 
                                                    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> &person_likes_comment_r,
                                                    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> &person_likes_post_r,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_hasCreator_person_r, 
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> &post_hasCreator_person_r,
                                                    std::unordered_map<std::string,std::string> &comment_hasCreator_person, 
                                                    std::unordered_map<std::string,std::string> &post_hasCreator_person,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> &post_containerOf_forum,
                                                    std::unordered_map<std::string,std::string> &comment_0_0,
                                                    std::unordered_map<std::string,std::string> &post_0_0,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_comment_r,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_post_r,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_comment,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_post,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_hasTag_tag,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> &post_hasTag_tag,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> &comment_isLocatedIn_place,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> &post_isLocatedIn_place) {
    // 构建带权图
    Graph g;
    createLGraph(person_0_0,comment_hasCreator_person_r,post_hasCreator_person_r,comment_hasCreator_person,post_hasCreator_person,comment_replyOf_comment,comment_replyOf_post,g);
    // printMemoryUsage();
    GetCurrentProcessMemoryUsage();
    auto start = std::chrono::high_resolution_clock::now();
    // std::unordered_map<std::string,std::unordered_set<std::string>> person_message_set = merge_maps(comment_hasCreator_person_r, post_hasCreator_person_r);
    std::unordered_set<std::string> communityMessageSet; //存储当前社区中person发布的message，用于判断阈值
    unordered_map<string, std::unordered_set<std::string>> communitiesSet;  //存储每个社区中person
    std::unordered_set<std::string> pPerson;  //存储已分配到社区中person
    int communityNo=0;
    int thresh = 0;
    // 分配社区中的person
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
        // communityMessageSet.insert(person_message_set[center].begin(),person_message_set[center].end());
        communityMessageSet.insert(comment_hasCreator_person_r[center].begin(),comment_hasCreator_person_r[center].end());
        communityMessageSet.insert(post_hasCreator_person_r[center].begin(),post_hasCreator_person_r[center].end());
        thresh = communityMessageSet.size();
        if(thresh>=(messageNum/part)) {
            std::cout<<"divide community phase1 finished,community size:"<<communitiesSet[to_string(communityNo)].size()<<std::endl;
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
            // communityMessageSet.insert(person_message_set[tmp].begin(),person_message_set[tmp].end());
            communityMessageSet.insert(comment_hasCreator_person_r[tmp].begin(),comment_hasCreator_person_r[tmp].end());
            communityMessageSet.insert(post_hasCreator_person_r[tmp].begin(),post_hasCreator_person_r[tmp].end());
            thresh = communityMessageSet.size();
            if(thresh>=(messageNum/part)) {
                std::cout<<"divide community phase1 finished,community size:"<<communitiesSet[to_string(communityNo)].size()<<std::endl;
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
            // communityMessageSet.insert(person_message_set[cur].begin(),person_message_set[cur].end());
            communityMessageSet.insert(comment_hasCreator_person_r[cur].begin(),comment_hasCreator_person_r[cur].end());
            communityMessageSet.insert(post_hasCreator_person_r[cur].begin(),post_hasCreator_person_r[cur].end());
            thresh = communityMessageSet.size();
            if(thresh>=(messageNum/part)) {
                std::cout<<"divide community phase2 finished,community size:"<<communitiesSet[to_string(communityNo)].size()<<std::endl;
                communityNo++;
                communityMessageSet.clear();
                thresh = 0;
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout<<"divide community finished,uses "<<duration.count()<<"s"<<endl;
    // 将message按照personid%part划分到part集合中
    std::vector<std::unordered_set<std::string>> commentSetArray(part);
    std::vector<std::unordered_set<std::string>> postSetArray(part);
    std::string outPath="./output";
    std::string outPath1="./route";
    if (!createDirectory(outPath)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    outPath="./output/messages";
    if (!createDirectory(outPath)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    if (!createDirectory(outPath1)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    std::ofstream message2partf("./route/message2part.csv");
    std::ofstream person2partf("./route/person2messagePart.csv");
    int commuCnt = 0;
    // int mc_num = part/communityNo;// 每个社区的message分片数量
    for(auto community:communitiesSet) {
        for(auto person:community.second) {
            int tmpCnt=0;
            // while(tmpCnt++<mc_num) person2partf<<person<<'|'<<std::to_string(commuCnt * mc_num + tmpCnt)<<std::endl;
            person2partf<<person<<'|'<<community.first<<std::endl;
            // 将person的hasCreator的中心message分片
            for(auto i: comment_hasCreator_person_r[person]) {
                // long long int offset = std::stol(i)%mc_num;
                // commentSetArray[commuCnt * mc_num+offset].insert(i);
                // message2partf<<i<<'|'<<commuCnt * mc_num+offset<<std::endl;
                commentSetArray[std::stol(community.first)].insert(i);
                message2partf<<i<<'|'<<community.first<<std::endl;
            }
            for(auto i: post_hasCreator_person_r[person]) {
                // long long int offset = std::stol(i)%mc_num;
                // postSetArray[commuCnt * mc_num+offset].insert(i);
                // message2partf<<i<<'|'<<commuCnt * mc_num+offset<<std::endl;
                postSetArray[std::stol(community.first)].insert(i);
                message2partf<<i<<'|'<<community.first<<std::endl;
            }
        }
        // commuCnt++;
    }
    message2partf.close();
    person2partf.close();
    
    for(int i=0;i<communitiesSet.size();i++){
        // GetCurrentProcessMemoryUsage();
        outPath="./output/messages/";
        outPath+=std::to_string(i);
        // 如果文件夹不存在，则创建它
        if (!createDirectory(outPath)) {
            std::cerr << "Failed or no need to create directory."<<outPath << std::endl;
        }
        std::cout<<outPath<<std::endl;
        std::ofstream outfile(outPath+"/comment_hasCreator_person_0_0.csv");
        std::ofstream outfile0(outPath+"/post_hasCreator_person_0_0.csv");
        std::ofstream outfile1(outPath+"/comment_0_0.csv");
        std::ofstream outfile1_1(outPath+"/comment_0_0_added.csv");
        std::ofstream outfile2(outPath+"/post_0_0.csv");
        std::ofstream outfile2_1(outPath+"/post_0_0_added.csv");

        std::ofstream file3_1(outPath+"/comment_replyOf_comment_0_0.csv");
        std::ofstream file3_2(outPath+"/comment_replyOf_post_0_0.csv");

        
        std::ofstream outfile4(outPath+"/comment_isLocatedIn_place_0_0.csv");
        std::ofstream outfile5(outPath+"/post_isLocatedIn_place_0_0.csv");
        std::ofstream outfile6(outPath+"/comment_hasTag_tag_0_0.csv");
        std::ofstream outfile7(outPath+"/post_hasTag_tag_0_0.csv");

        std::ofstream outfile8(outPath+"/forum_containerOf_post_0_0.csv");
        std::ofstream outfile9(outPath+"/forum_0_0.csv");

        std::ofstream outfile10(outPath+"/person_0_0.csv");

        std::ofstream likes1(outPath+"/person_likes_comment_0_0.csv");
        std::ofstream likes2(outPath+"/person_likes_post_0_0.csv");

        std::ofstream person_comment(outPath+"/person-comment_0_0.csv");
        std::ofstream person_post(outPath+"/person-post_0_0.csv");
        // std::cout<<"open"<<std::endl;
        // GetCurrentProcessMemoryUsage();

        if (outfile.is_open()) {
            std::cout << "文件成功打开，现在可以写入数据。" << std::endl;
            // 在这里写入数据到文件
        } else {
            std::cout << "文件打开失败，请检查路径和权限：" << outPath+"/post_hasCreator_person_0_0.csv" << std::endl;
        }
        // 写入person
        for(auto i:communitiesSet[std::to_string(i)]) {
            outfile10<<i<<std::endl;
        }
        outfile10.close();

        
        // 去重comment2comment
        std::unordered_map<std::string,std::unordered_map<std::string,int>> visited; 
        // 去重comment2post
        std::unordered_map<std::string,std::unordered_map<std::string,int>> visited1; 
        // 去重当前分片的comment、post、forum
        std::unordered_set<std::string> tmpComment;
        std::unordered_set<std::string> tmpPost;
        std::unordered_set<std::string> tmpForum;
        // 存储commentPath
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
                        file3_1<<kk<<'|'<<k<<std::endl;
                }
            }
            // comment的后向comment和边
            tmpComment.insert(comment_replyOf_comment[k].begin(),comment_replyOf_comment[k].end());
            for(auto kk:comment_replyOf_comment[k]){
                if(visited[k].find(kk)==visited[k].end()) {
                        visited[k][kk] = 1;
                        file3_1<<k<<'|'<<kk<<std::endl;
                }
            }
            // comment的后向post和边
            tmpPost.insert(comment_replyOf_post[k].begin(),comment_replyOf_post[k].end());
            for(auto kk:comment_replyOf_post[k]){
                if(visited1[k].find(kk)==visited1[k].end()) {
                        visited1[k][kk] = 1;
                        file3_2<<k<<'|'<<kk<<std::endl;
                }
            }
        }
        for(auto k:postSetArray[i]) {
            // post的前向
            tmpComment.insert(comment_replyOf_post_r[k].begin(),comment_replyOf_post_r[k].end());
            for(auto kk:comment_replyOf_post_r[k]){
                if(visited1[kk].find(k)==visited1[kk].end()) {
                        visited1[kk][k] = 1;
                        file3_2<<kk<<'|'<<k<<std::endl;
                }
            }
        }
        
        // std::unordered_map<std::string,std::vector<std::string>> commentPath;
        for (auto k:commentSetArray[i]) {
            std::vector<std::string> tmp;
            dfs(k, comment_replyOf_comment, tmp, commentPath);
        }
        for(auto k=commentPath.begin();k!=commentPath.end();k++){
            std::string from = k->first;
            std::vector<std::string> tmp = k->second;
            if(tmp.size()==1) {
                // comment没有comment邻居
                // ！！！这里也要存post
                tmpComment.insert(tmp[0]); // 里面本来就存了
                // file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                for(std::string j : comment_replyOf_post[tmp[0]]) {
                    tmpPost.insert(j);
                    if(visited1[tmp[0]].find(j)==visited1[tmp[0]].end()) {
                        visited1[tmp[0]][j] = 1;
                        file3_2<<tmp[0]<<'|'<<j<<std::endl;
                    }
                    // file3_2<<tmp[0]<<'|'<<j<<std::endl;
                }
            } else if(tmp.size()==2) {
                // comment只有一步comment
                tmpComment.insert(tmp[1]);
                if(visited[tmp[0]].find(tmp[1])==visited[tmp[0]].end()) {
                    visited[tmp[0]][tmp[1]] = 1;
                    file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                }
                // file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                for(std::string j : comment_replyOf_post[tmp[1]]) {
                    // comment存在comment链
                    tmpPost.insert(j);
                    if(visited1[tmp[1]].find(j)==visited1[tmp[1]].end()) {
                        visited1[tmp[1]][j] = 1;
                        file3_2<<tmp[1]<<'|'<<j<<std::endl;
                    }
                    // file3_2<<tmp[1]<<'|'<<j<<std::endl;
                }
            } else {
                // 顶点有多步邻居
                tmpComment.insert(tmp[1]);
                if(visited[tmp[0]].find(tmp[1])==visited[tmp[0]].end()) {
                    visited[tmp[0]][tmp[1]] = 1;
                    file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                }
                // 这里把末端邻居合并为post
                for(std::string j : comment_replyOf_post[tmp[2]]) {
                    tmpPost.insert(j);
                    if(visited1[tmp[1]].find(j)==visited1[tmp[1]].end()) {
                        visited1[tmp[1]][j] = 1;
                        file3_2<<tmp[1]<<'|'<<j<<std::endl;
                    }
                    // 压缩了路径，不过好像导致comment2post变多了？
                }
            }
        }
        
        // std::cout<<"message extend finished"<<std::endl;
        // GetCurrentProcessMemoryUsage();

        // 写入comment
        for(auto k:tmpComment){
            // 这里改成中心comment存整行，非中心comment存id（需要加空列）
            // outfile1原版本
            // outfile1<<comment_0_0[k]<<std::endl;
            // outfile1新版本
            if(commentSetArray[i].find(k)!=commentSetArray[i].end()) {
                outfile<<k<<'|'<<comment_hasCreator_person[k]<<std::endl;
                // 是中心message
                outfile1<<comment_0_0[k]<<std::endl;
                person_comment<<k<<std::endl;
                // likes
                for(auto j:person_likes_comment_r[k]){
                    likes1<<j.first<<'|'<<k<<'|'<<j.second<<std::endl;
                }
            } else {
                outfile1_1<<k<<std::endl;
            }
            for(auto j:comment_isLocatedIn_place[k]){
                outfile4<<k<<'|'<<j<<std::endl;
            }
            for(auto j:comment_hasTag_tag[k]){
                outfile6<<k<<'|'<<j<<std::endl;
            }
        }
        // 写入post
        for(auto k:tmpPost){
            // outfile2<<post_0_0[k]<<std::endl;
            if(postSetArray[i].find(k)!=postSetArray[i].end()) {
                outfile0<<k<<'|'<<post_hasCreator_person[k]<<std::endl;
                // 是中心message
                outfile2<<post_0_0[k]<<std::endl;
                person_post<<k<<std::endl;
                // likes
                for(auto j:person_likes_post_r[k]){
                    likes2<<j.first<<'|'<<k<<'|'<<j.second<<std::endl;
                }
            } else {
                outfile2_1<<k<<std::endl;
            }
            for(auto j:post_isLocatedIn_place[k]){
                outfile5<<k<<'|'<<j<<std::endl;
            }
            for(auto j:post_hasTag_tag[k]){
                outfile7<<k<<'|'<<j<<std::endl;
            }
            for(auto j:post_containerOf_forum[k]){
                outfile8<<j<<'|'<<k<<std::endl;
                tmpForum.insert(j);
            }
        }
        // std::cout<<"forum extend finished"<<std::endl;
        // GetCurrentProcessMemoryUsage();
        // 写入forum
        for(auto k:tmpForum){
            outfile9<<k<<std::endl;
        }
        

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
        
        outfile.close();
        outfile0.close();
        outfile1.close();
        outfile1_1.close();
        outfile2.close();
        outfile2_1.close();
        file3_1.close();
        file3_2.close();
        outfile4.close();
        outfile5.close();
        outfile6.close();
        outfile7.close();
        outfile8.close();
        outfile9.close();
        likes1.close();
        likes2.close();
        person_comment.close();
        person_post.close();
        // 强制释放内存
        // std::cout<<"close"<<std::endl;
        // GetCurrentProcessMemoryUsage();
        
    }
    // printMemoryUsage();
    communityMessageSet.clear();
    communitiesSet.clear();
    pPerson.clear();
}


// 多线程版本
void divideMessageCommunityThreads(int part, int messageNum, std::unordered_map<std::string,std::string> &person_0_0,
                                            std::unordered_map<std::string,std::unordered_map<std::string,std::string>> &person_likes_comment_r,
                                            std::unordered_map<std::string,std::unordered_map<std::string,std::string>> &person_likes_post_r,
                                            std::unordered_map<std::string,std::unordered_set<std::string>> &comment_hasCreator_person_r, 
                                            std::unordered_map<std::string,std::unordered_set<std::string>> &post_hasCreator_person_r,
                                            std::unordered_map<std::string,std::string> &comment_hasCreator_person, 
                                            std::unordered_map<std::string,std::string> &post_hasCreator_person,
                                            std::unordered_map<std::string,std::unordered_set<std::string>> &post_containerOf_forum,
                                            std::unordered_map<std::string,std::string> &comment_0_0,
                                            std::unordered_map<std::string,std::string> &post_0_0,
                                            std::unordered_map<std::string,std::string> &forum_0_0,
                                            std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_comment_r,
                                            std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_post_r,
                                            std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_comment,
                                            std::unordered_map<std::string,std::unordered_set<std::string>> &comment_replyOf_post,
                                            std::unordered_map<std::string,std::unordered_set<std::string>> &comment_hasTag_tag,
                                            std::unordered_map<std::string,std::unordered_set<std::string>> &post_hasTag_tag,
                                            std::unordered_map<std::string,std::unordered_set<std::string>> &comment_isLocatedIn_place,
                                            std::unordered_map<std::string,std::unordered_set<std::string>> &post_isLocatedIn_place) {
    auto start = std::chrono::high_resolution_clock::now();
    // 构建带权图
    Graph g;
    createLGraph(person_0_0,comment_hasCreator_person_r,post_hasCreator_person_r,comment_hasCreator_person,post_hasCreator_person,comment_replyOf_comment,comment_replyOf_post,g);
    auto LGraphEnd = std::chrono::high_resolution_clock::now();
    auto LGraphDuration = std::chrono::duration_cast<std::chrono::seconds>(LGraphEnd - start);
    std::cout << "person replyOf Graph created uses:" << LGraphDuration.count() << " seconds" << std::endl;
    std::unordered_map<std::string,std::unordered_set<std::string>> person_message_set = merge_maps(comment_hasCreator_person_r, post_hasCreator_person_r);
    unordered_map<string, std::unordered_set<std::string>> communityMessageSet; //存储社区中person发布的message
    unordered_map<string, std::unordered_set<std::string>> communitiesSet;  //存储社区中person
    std::unordered_set<std::string> pPerson;  //存储分配到社区中person
    int communityNo=0;
    int thresh = 0;
    // 分配社区中的person
    while(g.getAdjList().size()>0) {
        // 每次取第一个点作为社区中心
        string center = g.getAdjList().begin()->first;
        if(pPerson.find(center)!=pPerson.end()) {
            cout<<"same center"<<center<<endl;
        }
        pPerson.insert(center);
        thresh += person_message_set[center].size();
        communitiesSet[to_string(communityNo)].insert(center);
        communityMessageSet[to_string(communityNo)].insert(person_message_set[center].begin(),person_message_set[center].end());
        if(thresh>=(messageNum/part)) {
            std::cout<<"divide community phase1 finished,community size:"<<communitiesSet[to_string(communityNo)].size()<<std::endl;
            communityNo++;
            thresh = 0;
            g.removeNodes({center});
            continue;
        }
        // 取出社区中心的最大邻居，直到社区超过阈值或者没有邻居
        while(g.getNeighbor(center).size()!=0) {
            string tmp = g.getMaxNeighbor(center);
            if(pPerson.find(tmp)!=pPerson.end()) {
                cout<<"same tmp"<<tmp<<endl;
            }
            pPerson.insert(tmp);
            thresh += person_message_set[tmp].size();
            communitiesSet[to_string(communityNo)].insert(tmp);
            communityMessageSet[to_string(communityNo)].insert(person_message_set[tmp].begin(),person_message_set[tmp].end());
            g.removeNodes({tmp});
            if(thresh>=(messageNum/part)) {
                std::cout<<"divide community phase1 finished,community size:"<<communitiesSet[to_string(communityNo)].size()<<std::endl;
                communityNo++;
                thresh = 0;
                break;
            }
        }
        g.removeNodes({center});
    }
    // 分配不在社区中person，就是没有发布message的person或者message没有关联的person
    for(auto i:person_0_0) {
        if(pPerson.find(i.first)==pPerson.end()) {
            string cur = i.first;
            pPerson.insert(cur);
            thresh += person_message_set[cur].size();
            communitiesSet[to_string(communityNo)].insert(cur);
            communityMessageSet[to_string(communityNo)].insert(person_message_set[cur].begin(),person_message_set[cur].end());
            if(thresh>=(messageNum/part)) {
                std::cout<<"divide community phase2 finished,community size:"<<communitiesSet[to_string(communityNo)].size()<<std::endl;
                communityNo++;
                thresh = 0;
            }
        }
    }
    auto CommunityEnd = std::chrono::high_resolution_clock::now();
    auto CommunityDuration = std::chrono::duration_cast<std::chrono::seconds>(CommunityEnd - LGraphEnd);
    std::cout << "Community divided uses:" << CommunityDuration.count() << " seconds" << std::endl;
    
    // 将message按照personid%part划分到part集合中
    std::vector<std::unordered_set<std::string>> commentSetArray(part);
    std::vector<std::unordered_set<std::string>> postSetArray(part);
    std::string outPath="./output";
    std::string outPath1="./route";
    if (!createDirectory(outPath)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    outPath="./output/messages";
    if (!createDirectory(outPath)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    if (!createDirectory(outPath1)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }

    SafeQueue message2partfqueue;
    std::thread message2partfThread(writer, std::ref(message2partfqueue), "./route/message2part.csv");

    SafeQueue person2partfqueue;
    std::thread person2partfThread(writer, std::ref(person2partfqueue), "./route/person2messagePart.csv");

    // std::ofstream message2partf("./route/message2part.csv");
    // std::ofstream person2partf("./route/person2messagePart.csv");
    
    int commuCnt = 0;
    // int mc_num = part/communityNo;// 每个社区的message分片数量
    for(auto community:communitiesSet) {
        for(auto person:community.second) {
            int tmpCnt=0;
            string person2partf=person+'|'+community.first+'\n';
            person2partfqueue.push(person2partf);
            // 将person的hasCreator的中心message分片
            for(auto i: comment_hasCreator_person_r[person]) {
                
                commentSetArray[std::stol(community.first)].insert(i);
                string message2partf=i+'|'+community.first+'\n';
                message2partfqueue.push(message2partf);
            }
            for(auto i: post_hasCreator_person_r[person]) {
                
                postSetArray[std::stol(community.first)].insert(i);
                string message2partf=i+'|'+community.first+'\n';
                message2partfqueue.push(message2partf);
            }
        }
        // commuCnt++;
    }
    person2partfqueue.setDone();
    message2partfqueue.setDone();
    person2partfThread.join();
    message2partfThread.join();
    
    for(int i=0;i<communitiesSet.size();i++){
        outPath="./output/messages/";
        outPath+=std::to_string(i);
        // 如果文件夹不存在，则创建它
        if (!createDirectory(outPath)) {
            std::cerr << "Failed or no need to create directory."<<outPath << std::endl;
        }
        std::cout<<outPath<<std::endl;
        // std::ofstream outfile(outPath+"/comment_hasCreator_person_0_0.csv");
        SafeQueue outfileQueue;
        std::thread outfileThread(writer, std::ref(outfileQueue), outPath+"/comment_hasCreator_person_0_0.csv");
        // std::ofstream outfile0(outPath+"/post_hasCreator_person_0_0.csv");
        SafeQueue outfile0Queue;
        std::thread outfile0Thread(writer, std::ref(outfile0Queue), outPath+"/post_hasCreator_person_0_0.csv");
        // std::ofstream outfile1(outPath+"/comment_0_0.csv");
        SafeQueue outfile1Queue;
        std::thread outfile1Thread(writer, std::ref(outfile1Queue), outPath+"/comment_0_0.csv");
        // std::ofstream outfile1_1(outPath+"/comment_0_0_added.csv");
        SafeQueue outfile11Queue;
        std::thread outfile11Thread(writer, std::ref(outfile11Queue), outPath+"/comment_0_0_added.csv");
        // std::ofstream outfile2(outPath+"/post_0_0.csv");
        SafeQueue outfile2Queue;
        std::thread outfile2Thread(writer, std::ref(outfile2Queue), outPath+"/post_0_0.csv");
        // std::ofstream outfile2_1(outPath+"/post_0_0_added.csv");
        SafeQueue outfile21Queue;
        std::thread outfile21Thread(writer, std::ref(outfile21Queue), outPath+"/post_0_0_added.csv");

        // std::ofstream file3_1(outPath+"/comment_replyOf_comment_0_0.csv");
        SafeQueue outfile31Queue;
        std::thread outfile31Thread(writer, std::ref(outfile31Queue), outPath+"/comment_replyOf_comment_0_0.csv");
        // std::ofstream file3_2(outPath+"/comment_replyOf_post_0_0.csv");
        SafeQueue outfile32Queue;
        std::thread outfile32Thread(writer, std::ref(outfile32Queue), outPath+"/comment_replyOf_post_0_0.csv");

        
        // std::ofstream outfile4(outPath+"/comment_isLocatedIn_place_0_0.csv");
        SafeQueue outfile4Queue;
        std::thread outfile4Thread(writer, std::ref(outfile4Queue), outPath+"/comment_isLocatedIn_place_0_0.csv");
        // std::ofstream outfile5(outPath+"/post_isLocatedIn_place_0_0.csv");
        SafeQueue outfile5Queue;
        std::thread outfile5Thread(writer, std::ref(outfile5Queue), outPath+"/post_isLocatedIn_place_0_0.csv");
        // std::ofstream outfile6(outPath+"/comment_hasTag_tag_0_0.csv");
        SafeQueue outfile6Queue;
        std::thread outfile6Thread(writer, std::ref(outfile6Queue), outPath+"/comment_hasTag_tag_0_0.csv");
        // std::ofstream outfile7(outPath+"/post_hasTag_tag_0_0.csv");
        SafeQueue outfile7Queue;
        std::thread outfile7Thread(writer, std::ref(outfile7Queue), outPath+"/post_hasTag_tag_0_0.csv");

        // std::ofstream outfile8(outPath+"/forum_containerOf_post_0_0.csv");
        SafeQueue outfile8Queue;
        std::thread outfile8Thread(writer, std::ref(outfile8Queue), outPath+"/forum_containerOf_post_0_0.csv");
        // std::ofstream outfile9(outPath+"/forum_0_0.csv");
        SafeQueue outfile9Queue;
        std::thread outfile9Thread(writer, std::ref(outfile9Queue), outPath+"/forum_0_0.csv");

        // std::ofstream outfile10(outPath+"/person_0_0.csv");
        SafeQueue outfile10Queue;
        std::thread outfile10Thread(writer, std::ref(outfile10Queue), outPath+"/person_0_0.csv");

        // std::ofstream outfile10(outPath+"/person_0_0.csv");
        SafeQueue likes1Queue;
        std::thread likes1Thread(writer, std::ref(likes1Queue), outPath+"/person_likes_comment_0_0.csv");

        // std::ofstream outfile10(outPath+"/person_0_0.csv");
        SafeQueue likes2Queue;
        std::thread likes2Thread(writer, std::ref(likes2Queue), outPath+"/person_likes_post_0_0.csv");

        // if (outfile.is_open()) {
        //     std::cout << "文件成功打开，现在可以写入数据。" << std::endl;
        //     // 在这里写入数据到文件
        // } else {
        //     std::cout << "文件打开失败，请检查路径和权限：" << outPath+"/post_hasCreator_person_0_0.csv" << std::endl;
        // }
        // 写入person
        for(auto i:communitiesSet[std::to_string(i)]) {
            // outfile10<<i<<std::endl;
            string input = i+'\n';
            outfile10Queue.push(input);
        }
        outfile10Queue.setDone();
        // outfile10.close();

        
        // 去重comment2comment
        std::unordered_map<std::string,std::unordered_map<std::string,int>> visited; 
        // 去重comment2post
        std::unordered_map<std::string,std::unordered_map<std::string,int>> visited1; 
        
        std::set<std::string> tmpComment;
        std::set<std::string> tmpPost;
        std::set<std::string> tmpForum;
        std::set<std::string> tmpPerson;
        tmpComment.insert(commentSetArray[i].begin(),commentSetArray[i].end());
        tmpPost.insert(postSetArray[i].begin(),postSetArray[i].end());

        // 开始划分message
        for(auto k:commentSetArray[i]) {
            // comment的前向comment和边
            tmpComment.insert(comment_replyOf_comment_r[k].begin(),comment_replyOf_comment_r[k].end());
            for(auto kk:comment_replyOf_comment_r[k]){
                if(visited[kk].find(k)==visited[kk].end()) {
                        visited[kk][k] = 1;
                        // file3_1<<kk<<'|'<<k<<std::endl;
                        string input = kk+'|'+k+'\n';
                        outfile31Queue.push(input);
                }
            }
            // comment的后向comment和边
            tmpComment.insert(comment_replyOf_comment[k].begin(),comment_replyOf_comment[k].end());
            for(auto kk:comment_replyOf_comment[k]){
                if(visited[k].find(kk)==visited[k].end()) {
                        visited[k][kk] = 1;
                        // file3_1<<k<<'|'<<kk<<std::endl;
                        string input = k+'|'+kk+'\n';
                        outfile31Queue.push(input);
                }
            }
            // comment的后向post和边
            tmpPost.insert(comment_replyOf_post[k].begin(),comment_replyOf_post[k].end());
            for(auto kk:comment_replyOf_post[k]){
                if(visited1[k].find(kk)==visited1[k].end()) {
                        visited1[k][kk] = 1;
                        // file3_2<<k<<'|'<<kk<<std::endl;
                        string input = k+'|'+kk+'\n';
                        outfile32Queue.push(input);
                }
            }
        }
        for(auto k:postSetArray[i]) {
            // post的前向
            tmpComment.insert(comment_replyOf_post_r[k].begin(),comment_replyOf_post_r[k].end());
            for(auto kk:comment_replyOf_post_r[k]){
                if(visited1[kk].find(k)==visited1[kk].end()) {
                        visited1[kk][k] = 1;
                        // file3_2<<kk<<'|'<<k<<std::endl;
                        string input = kk+'|'+k+'\n';
                        outfile32Queue.push(input);
                }
            }
        }
        
        std::unordered_map<std::string,std::vector<std::string>> commentPath;
        for (auto k:commentSetArray[i]) {
            std::vector<std::string> tmp;
            dfs(k, comment_replyOf_comment, tmp, commentPath);
        }
        for(auto k=commentPath.begin();k!=commentPath.end();k++){
            std::string from = k->first;
            std::vector<std::string> tmp = k->second;
            if(tmp.size()==1) {
                // comment没有comment邻居
                // ！！！这里也要存post
                tmpComment.insert(tmp[0]); // 里面本来就存了
                for(std::string j : comment_replyOf_post[tmp[0]]) {
                    tmpPost.insert(j);
                    if(visited1[tmp[0]].find(j)==visited1[tmp[0]].end()) {
                        visited1[tmp[0]][j] = 1;
                        // file3_2<<tmp[0]<<'|'<<j<<std::endl;
                        string input = tmp[0]+'|'+j+'\n';
                        outfile32Queue.push(input);
                    }
                }
            } else if(tmp.size()==2) {
                // comment只有一步comment
                tmpComment.insert(tmp[1]);
                if(visited[tmp[0]].find(tmp[1])==visited[tmp[0]].end()) {
                    visited[tmp[0]][tmp[1]] = 1;
                    // file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                    string input = tmp[0]+'|'+tmp[1]+'\n';
                    outfile31Queue.push(input);
                }
                for(std::string j : comment_replyOf_post[tmp[1]]) {
                    // comment存在comment链
                    tmpPost.insert(j);
                    if(visited1[tmp[1]].find(j)==visited1[tmp[1]].end()) {
                        visited1[tmp[1]][j] = 1;
                        // file3_2<<tmp[1]<<'|'<<j<<std::endl;
                        string input = tmp[1]+'|'+j+'\n';
                        outfile32Queue.push(input);
                    }
                }
            } else {
                // 顶点有多步邻居
                tmpComment.insert(tmp[1]);
                if(visited[tmp[0]].find(tmp[1])==visited[tmp[0]].end()) {
                    visited[tmp[0]][tmp[1]] = 1;
                    // file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                    string input = tmp[0]+'|'+tmp[1]+'\n';
                    outfile31Queue.push(input);
                }
                // 这里把末端邻居合并为post
                for(std::string j : comment_replyOf_post[tmp[2]]) {
                    tmpPost.insert(j);
                    if(visited1[tmp[1]].find(j)==visited1[tmp[1]].end()) {
                        visited1[tmp[1]][j] = 1;
                        // file3_2<<tmp[1]<<'|'<<j<<std::endl;
                        string input = tmp[1]+'|'+j+'\n';
                        outfile32Queue.push(input);
                    }
                    // 压缩了路径，不过好像导致comment2post变多了？
                }
            }
        }
        outfile31Queue.setDone();
        outfile32Queue.setDone();

        // 写入comment
        for(auto k:tmpComment){
            // 这里改成中心comment存整行，非中心comment存id（需要加空列）
            // outfile1原版本
            // outfile1<<comment_0_0[k]<<std::endl;
            // outfile1新版本
            if(commentSetArray[i].find(k)!=commentSetArray[i].end()) {
                // outfile<<k<<'|'<<comment_hasCreator_person[k]<<std::endl;
                string input = k+'|'+comment_hasCreator_person[k]+'\n';
                outfileQueue.push(input);
                // 是中心message
                // outfile1<<comment_0_0[k]<<std::endl;
                string input1 = comment_0_0[k]+'\n';
                outfile1Queue.push(input1);
                // likes
                for(auto j:person_likes_comment_r[k]){
                    // likes1<<j.first<<'|'<<k<<'|'<<j.second<<std::endl;
                    string inputLikes = j.first+'|'+k+'|'+j.second+'\n';
                    likes1Queue.push(inputLikes);
                }
            } else {
                // outfile1_1<<k<<std::endl;
                string input = k+'\n';
                outfile11Queue.push(input);
            }
            for(auto j:comment_isLocatedIn_place[k]){
                // outfile4<<k<<'|'<<j<<std::endl;
                string input = k+'|'+j+'\n';
                outfile4Queue.push(input);
            }
            for(auto j:comment_hasTag_tag[k]){
                // outfile6<<k<<'|'<<j<<std::endl;
                string input = k+'|'+j+'\n';
                outfile6Queue.push(input);
            }
        }
        likes1Queue.setDone();
        outfileQueue.setDone();
        outfile1Queue.setDone();
        outfile11Queue.setDone();
        outfile4Queue.setDone();
        outfile6Queue.setDone();
        // 写入post
        for(auto k:tmpPost){
            // outfile2<<post_0_0[k]<<std::endl;
            if(postSetArray[i].find(k)!=postSetArray[i].end()) {
                // outfile0<<k<<'|'<<post_hasCreator_person[k]<<std::endl;
                string input = k+'|'+post_hasCreator_person[k]+'\n';
                outfile0Queue.push(input);
                // 是中心message
                // outfile2<<post_0_0[k]<<std::endl;
                string input1 = post_0_0[k]+'\n';
                outfile2Queue.push(input1);
                // likes
                for(auto j:person_likes_post_r[k]){
                    // likes2<<j.first<<'|'<<k<<'|'<<j.second<<std::endl;
                    string inputLikes = j.first+'|'+k+'|'+j.second+'\n';
                    likes2Queue.push(inputLikes);
                }
            } else {
                // outfile2_1<<k<<std::endl;
                string input = k+'\n';
                outfile21Queue.push(input);
            }
            for(auto j:post_isLocatedIn_place[k]){
                // outfile5<<k<<'|'<<j<<std::endl;
                string input = k+'|'+j+'\n';
                outfile5Queue.push(input);
            }
            for(auto j:post_hasTag_tag[k]){
                // outfile7<<k<<'|'<<j<<std::endl;
                string input = k+'|'+j+'\n';
                outfile7Queue.push(input);
            }
            for(auto j:post_containerOf_forum[k]){
                // outfile8<<j<<'|'<<k<<std::endl;
                string input = j+'|'+k+'\n';
                outfile8Queue.push(input);
                tmpForum.insert(j);
            }
        }
        likes2Queue.setDone();
        outfile0Queue.setDone();
        outfile2Queue.setDone();
        outfile21Queue.setDone();
        outfile5Queue.setDone();
        outfile7Queue.setDone();
        outfile8Queue.setDone();
        // 写入forum
        for(auto k:tmpForum){
            // outfile9<<k<<std::endl;
            string input = k+'\n';
            outfile9Queue.push(input);
        }
        outfile9Queue.setDone();
        
        // outfile.close();
        outfileThread.join();
        // outfile0.close();
        outfile0Thread.join();
        // outfile1.close();
        outfile1Thread.join();
        // outfile1_1.close();
        outfile11Thread.join();
        // outfile2.close();
        outfile2Thread.join();
        // outfile2_1.close();
        outfile21Thread.join();
        // file3_1.close();
        outfile31Thread.join();
        // file3_2.close();
        outfile32Thread.join();
        // outfile4.close();
        outfile4Thread.join();
        // outfile5.close();
        outfile5Thread.join();
        // outfile6.close();
        outfile6Thread.join();
        // outfile7.close();
        outfile7Thread.join();
        // outfile8.close();
        outfile8Thread.join();
        // outfile9.close();
        outfile9Thread.join();
        outfile10Thread.join();

        likes1Thread.join();
        likes2Thread.join();
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - CommunityEnd);
    std::cout << "Message divided uses:" << duration.count() << " seconds" << std::endl;

}

void divideMessageHashAll(int part, std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasCreator_person_r, 
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasCreator_person_r,
                                                    std::unordered_map<std::string,std::string> comment_hasCreator_person, 
                                                    std::unordered_map<std::string,std::string> post_hasCreator_person,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> post_containerOf_forum,
                                                    std::unordered_map<std::string,std::string> comment_0_0,
                                                    std::unordered_map<std::string,std::string> post_0_0,
                                                    std::unordered_map<std::string,std::string> forum_0_0,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> comment_replyOf_comment_r,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> comment_replyOf_post_r,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> comment_replyOf_comment,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> comment_replyOf_post,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> comment_hasTag_tag,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> post_hasTag_tag,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> comment_isLocatedIn_place,
                                                    std::unordered_map<std::string,std::unordered_set<std::string>> post_isLocatedIn_place) {
    // 将message按照personid%part划分到part集合中
    std::vector<std::unordered_set<std::string>> commentSetArray(part);
    std::vector<std::unordered_set<std::string>> postSetArray(part);
    std::string outPath="./output";
    std::string outPath1="./route";
    if (!createDirectory(outPath)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    outPath="./output/messages";
    if (!createDirectory(outPath)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    if (!createDirectory(outPath1)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    std::ofstream message2partf("./route/message2part.csv");
    // 对messageid哈希,按照messageid%part划分，需要添加路由表message2partf
    for(auto i: comment_0_0) {
        long long int index = std::stol(i.first)%part;
        commentSetArray[index].insert(i.first);
        message2partf<<i.first<<'|'<<index<<std::endl;
    }
    for(auto i: post_0_0) {
        long long int index = std::stol(i.first)%part;
        postSetArray[index].insert(i.first);
        message2partf<<i.first<<'|'<<index<<std::endl;
    }
    message2partf.close();
    std::unordered_map<std::string,std::unordered_set<std::string>> person2messagePart;
    for(int i=0;i<part;i++){
        outPath="./output/messages/";
        outPath+=std::to_string(i);
        // 如果文件夹不存在，则创建它
        if (!createDirectory(outPath)) {
            std::cerr << "Failed or no need to create directory."<<outPath << std::endl;
        }
        std::cout<<outPath<<std::endl;
        std::ofstream outfile(outPath+"/comment_hasCreator_person_0_0.csv");
        std::ofstream outfile0(outPath+"/post_hasCreator_person_0_0.csv");
        std::ofstream outfile1(outPath+"/comment_0_0.csv");
        std::ofstream outfile1_1(outPath+"/comment_0_0_added.csv");
        std::ofstream outfile2(outPath+"/post_0_0.csv");
        std::ofstream outfile2_1(outPath+"/post_0_0_added.csv");

        std::ofstream file3_1(outPath+"/comment_replyOf_comment_0_0.csv");
        std::ofstream file3_2(outPath+"/comment_replyOf_post_0_0.csv");

        
        std::ofstream outfile4(outPath+"/comment_isLocatedIn_place_0_0.csv");
        std::ofstream outfile5(outPath+"/post_isLocatedIn_place_0_0.csv");
        std::ofstream outfile6(outPath+"/comment_hasTag_tag_0_0.csv");
        std::ofstream outfile7(outPath+"/post_hasTag_tag_0_0.csv");

        std::ofstream outfile8(outPath+"/forum_containerOf_post_0_0.csv");
        std::ofstream outfile9(outPath+"/forum_0_0.csv");

        if (outfile.is_open()) {
            std::cout << "文件成功打开，现在可以写入数据。" << std::endl;
            // 在这里写入数据到文件
        } else {
            std::cout << "文件打开失败，请检查路径和权限：" << outPath+"/post_hasCreator_person_0_0.csv" << std::endl;
        }
        // 去重comment2comment
        std::unordered_map<std::string,std::unordered_map<std::string,int>> visited; 
        // 去重comment2post
        std::unordered_map<std::string,std::unordered_map<std::string,int>> visited1; 
        
        std::unordered_set<std::string> tmpComment;
        std::unordered_set<std::string> tmpPost;
        std::unordered_set<std::string> tmpForum;
        tmpComment.insert(commentSetArray[i].begin(),commentSetArray[i].end());
        tmpPost.insert(postSetArray[i].begin(),postSetArray[i].end());

        // 将person的hasCreator的中心message分片插入person2messagePart
        for(auto k:commentSetArray[i]) {
            person2messagePart[comment_hasCreator_person[k]].insert(std::to_string(i));
        }
        for(auto k:postSetArray[i]) {
            person2messagePart[post_hasCreator_person[k]].insert(std::to_string(i));
        }
        // 开始划分message
        for(auto k:commentSetArray[i]) {
            // comment的前向comment和边
            tmpComment.insert(comment_replyOf_comment_r[k].begin(),comment_replyOf_comment_r[k].end());
            for(auto kk:comment_replyOf_comment_r[k]){
                if(visited[kk].find(k)==visited[kk].end()) {
                        visited[kk][k] = 1;
                        file3_1<<kk<<'|'<<k<<std::endl;
                }
            }
            // comment的后向comment和边
            tmpComment.insert(comment_replyOf_comment[k].begin(),comment_replyOf_comment[k].end());
            for(auto kk:comment_replyOf_comment[k]){
                if(visited[k].find(kk)==visited[k].end()) {
                        visited[k][kk] = 1;
                        file3_1<<k<<'|'<<kk<<std::endl;
                }
            }
            // comment的后向post和边
            tmpPost.insert(comment_replyOf_post[k].begin(),comment_replyOf_post[k].end());
            for(auto kk:comment_replyOf_post[k]){
                if(visited1[k].find(kk)==visited1[k].end()) {
                        visited1[k][kk] = 1;
                        file3_2<<k<<'|'<<kk<<std::endl;
                }
            }
        }
        for(auto k:postSetArray[i]) {
            // post的前向
            tmpComment.insert(comment_replyOf_post_r[k].begin(),comment_replyOf_post_r[k].end());
            for(auto kk:comment_replyOf_post_r[k]){
                if(visited1[kk].find(k)==visited1[kk].end()) {
                        visited1[kk][k] = 1;
                        file3_2<<kk<<'|'<<k<<std::endl;
                }
            }
        }
        
        std::unordered_map<std::string,std::vector<std::string>> commentPath;
        for (auto k:commentSetArray[i]) {
            std::vector<std::string> tmp;
            dfs(k, comment_replyOf_comment, tmp, commentPath);
        }
        for(auto k=commentPath.begin();k!=commentPath.end();k++){
            std::string from = k->first;
            std::vector<std::string> tmp = k->second;
            if(tmp.size()==1) {
                // comment没有comment邻居
                // ！！！这里也要存post
                tmpComment.insert(tmp[0]); // 里面本来就存了
                // file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                for(std::string j : comment_replyOf_post[tmp[0]]) {
                    tmpPost.insert(j);
                    if(visited1[tmp[0]].find(j)==visited1[tmp[0]].end()) {
                        visited1[tmp[0]][j] = 1;
                        file3_2<<tmp[0]<<'|'<<j<<std::endl;
                    }
                    // file3_2<<tmp[0]<<'|'<<j<<std::endl;
                }
            } else if(tmp.size()==2) {
                // comment只有一步comment
                tmpComment.insert(tmp[1]);
                if(visited[tmp[0]].find(tmp[1])==visited[tmp[0]].end()) {
                    visited[tmp[0]][tmp[1]] = 1;
                    file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                }
                // file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                for(std::string j : comment_replyOf_post[tmp[1]]) {
                    // comment存在comment链
                    tmpPost.insert(j);
                    if(visited1[tmp[1]].find(j)==visited1[tmp[1]].end()) {
                        visited1[tmp[1]][j] = 1;
                        file3_2<<tmp[1]<<'|'<<j<<std::endl;
                    }
                    // file3_2<<tmp[1]<<'|'<<j<<std::endl;
                }
            } else {
                // 顶点有多步邻居
                tmpComment.insert(tmp[1]);
                if(visited[tmp[0]].find(tmp[1])==visited[tmp[0]].end()) {
                    visited[tmp[0]][tmp[1]] = 1;
                    file3_1<<tmp[0]<<'|'<<tmp[1]<<std::endl;
                }
                // 这里把末端邻居合并为post
                for(std::string j : comment_replyOf_post[tmp[2]]) {
                    tmpPost.insert(j);
                    if(visited1[tmp[1]].find(j)==visited1[tmp[1]].end()) {
                        visited1[tmp[1]][j] = 1;
                        file3_2<<tmp[1]<<'|'<<j<<std::endl;
                    }
                    // 压缩了路径，不过好像导致comment2post变多了？
                }
            }
        }

        // 写入comment
        for(auto k:tmpComment){
            // 这里改成中心comment存整行，非中心comment存id（需要加空列）
            // outfile1原版本
            // outfile1<<comment_0_0[k]<<std::endl;
            // outfile1新版本
            if(commentSetArray[i].find(k)==commentSetArray[i].end()) {
                outfile<<k<<'|'<<comment_hasCreator_person[k]<<std::endl;
                // 是中心message
                outfile1<<comment_0_0[k]<<std::endl;
            } else {
                outfile1_1<<k<<std::endl;
            }
            for(auto j:comment_isLocatedIn_place[k]){
                outfile4<<k<<'|'<<j<<std::endl;
            }
            for(auto j:comment_hasTag_tag[k]){
                outfile6<<k<<'|'<<j<<std::endl;
            }
        }
        // 写入post
        for(auto k:tmpPost){
            // outfile2<<post_0_0[k]<<std::endl;
            if(postSetArray[i].find(k)==postSetArray[i].end()) {
                outfile0<<k<<'|'<<post_hasCreator_person[k]<<std::endl;
                // 是中心message
                outfile2<<post_0_0[k]<<std::endl;
            } else {
                outfile2_1<<k<<std::endl;
            }
            for(auto j:post_isLocatedIn_place[k]){
                outfile5<<k<<'|'<<j<<std::endl;
            }
            for(auto j:post_hasTag_tag[k]){
                outfile7<<k<<'|'<<j<<std::endl;
            }
            for(auto j:post_containerOf_forum[k]){
                outfile8<<j<<'|'<<k<<std::endl;
                tmpForum.insert(j);
            }
        }
        // 写入forum
        for(auto k:tmpForum){
            outfile9<<k<<std::endl;
        }
        
        outfile.close();
        outfile0.close();
        outfile1.close();
        outfile1_1.close();
        outfile2.close();
        outfile2_1.close();
        file3_1.close();
        file3_2.close();
        outfile4.close();
        outfile5.close();
        outfile6.close();
        outfile7.close();

        outfile8.close();
        outfile9.close();
    }

    std::ofstream outfile8("./route/person2messagePart.csv");
    for(auto k:person2messagePart){
        for(auto j:k.second){
            outfile8<<k.first<<'|'<<j<<std::endl;
        }
    }
    outfile8.close();

}

void onecol(std::unordered_map<std::string,std::unordered_map<std::string,std::string>> &forum_hasMember_person, std::unordered_map<std::string,std::unordered_set<std::string>> &forum_containerOf_post_r) {
    std::string outPath="./output/";
    std::ofstream outfile(outPath+"/forum_0_0-person.csv");
    std::ofstream outfile1(outPath+"/post_0_0-forum.csv");
    std::set<std::string> forum;
    std::set<std::string> post;
    for(auto i:forum_hasMember_person) {
        forum.insert(i.first);
    }
    for(auto i:forum_containerOf_post_r) {
        post.insert(i.first);
    }
    for(auto i:forum) {
        outfile<<i<<std::endl;
    }
    for(auto i:post) {
        outfile1<<i<<std::endl;
    }
}
void personAdd(//std::unordered_map<std::string,std::unordered_map<std::string,std::string>> &forum_hasMember_person,
                std::unordered_map<std::string, std::string> &comment_0_0,
                std::unordered_map<std::string, std::string> &post_0_0) {
    std::string outPath="./output/person";
    if (!createDirectory(outPath)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    // 移到forum分片
    // std::ofstream outfile(outPath+"/forum_0_0.csv");
    // std::set<std::string> forum;
    // for(auto i:forum_hasMember_person) {
    //     forum.insert(i.first);
    // }
    // for(auto i:forum) {
    //     outfile<<i<<std::endl;
    // }
    // outfile.close();

    std::ofstream outfile(outPath+"/comment_0_0.csv");
    for(auto i:comment_0_0) {
        outfile<<i.first<<std::endl;
    }
    outfile.close();
    
    outfile.open(outPath+"/post_0_0.csv");
    for(auto i:post_0_0) {
        outfile<<i.first<<std::endl;
    }
    outfile.close();
}


void forumAdd(std::unordered_map<std::string, std::string> &post_0_0) {
    std::string outPath="./output/forum";
    if (!createDirectory(outPath)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    std::ofstream outfile1(outPath+"/post_0_0.csv");
    for(auto i:post_0_0) {
        outfile1<<i.first<<std::endl;
    }
    outfile1.close();
}


void forumAdd1(std::unordered_map<std::string,std::unordered_set<std::string>> &forum_containerOf_post_r, std::unordered_map<std::string,std::string> &post_hasCreator_person) {
    std::string outPath="./output/forum";
    if (!createDirectory(outPath)) {
        // std::cerr << "Failed or no need to create directory." << std::endl;
    }
    std::ofstream outfile1(outPath+"/post_0_0.csv");
    // 这里forum是一个分片，其实就是全图
    std::ofstream outfile2(outPath+"/post_hasCreator_person_0_0.csv");
    std::set<std::string> post;
    std::set<std::string> person;
    // 注意，这里直接把反向边的post加入就不遍历forum了
    for(auto i:forum_containerOf_post_r) {
        post.insert(i.first);
    }
    for(auto i:post) {
        outfile1<<i<<std::endl;
        // post hasCreator person
        outfile2<<i<<'|'<<post_hasCreator_person[i]<<std::endl;
    }
    outfile1.close();
    outfile2.close();
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


int main() {
    // 获取开始时间点
    auto start = std::chrono::high_resolution_clock::now();

    // std::string initialDir = "/home/shared/data/social_network-csv_composite-longdateformatter-sf0.1/dynamic/";
    std::string initialDir = "/home/shared/data/social_network-csv_composite-longdateformatter-sf30/dynamic/";

    // 调用loadAll函数，加载所有数据
    std::unordered_map<std::string,std::unordered_set<std::string>> person_hasInterest_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> person_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_studyAt_organisation;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_workAt_organisation;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_likes_comment_r;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_likes_post_r;
    // std::unordered_map<std::string,std::unordered_map<std::string,std::string>> forum_hasMember_person_r;
    // std::unordered_map<std::string,std::unordered_map<std::string,std::string>> forum_hasMember_person;
    // std::unordered_map<std::string,std::unordered_set<std::string>> forum_hasModerator_person_r;
    // std::unordered_map<std::string,std::string> forum_hasModerator_person;
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
    // std::unordered_map<std::string, std::unordered_set<std::string>> reply1;
    // std::unordered_map<std::string, std::unordered_set<std::string>> reply2;
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

    
    int part = 100;
    Graph g;    //存储person-replyOf带权图
    unordered_map<string, std::unordered_set<std::string>> communitiesSet;
    std::vector<std::unordered_set<std::string>> commentSetArray(part);
    std::vector<std::unordered_set<std::string>> postSetArray(part);
    std::cout<<"initial memory:"<<std::endl;
    GetCurrentProcessMemoryUsage();

    // version1:
    // // origin loadAll
    // // loadAll(initialDir, person, person_hasInterest_tag, person_isLocatedIn_place, person_studyAt_organisation, person_workAt_organisation, person_likes_comment_r, person_likes_post_r, forum_hasMember_person_r, forum_hasMember_person, forum_hasModerator_person_r, forum_hasModerator_person, comment_hasCreator_person_r, comment_hasCreator_person, post_hasCreator_person_r, post_hasCreator_person, post_hasTag_tag, post_isLocatedIn_place, comment_hasTag_tag, comment_isLocatedIn_place, forum_containerOf_post, forum_containerOf_post_r, comment_replyOf_comment_r,comment_replyOf_post_r,comment_replyOf_comment,comment_replyOf_post, person_0_0, comment_0_0, post_0_0, forum_0_0);
    // new loadAll, remove useless load
    // loadAll(initialDir, person_hasInterest_tag, person_isLocatedIn_place, person_studyAt_organisation, person_workAt_organisation, person_likes_comment_r, person_likes_post_r, comment_hasCreator_person_r, comment_hasCreator_person, post_hasCreator_person_r, post_hasCreator_person, post_hasTag_tag, post_isLocatedIn_place, comment_hasTag_tag, comment_isLocatedIn_place, forum_containerOf_post, forum_containerOf_post_r, comment_replyOf_comment_r,comment_replyOf_post_r,comment_replyOf_comment,comment_replyOf_post, person_0_0, comment_0_0, post_0_0, forum_0_0);
    // auto loadEnd = std::chrono::high_resolution_clock::now();
    // auto loadDuration = std::chrono::duration_cast<std::chrono::seconds>(loadEnd - start);
    // std::cout << "load uses: " << loadDuration.count() << " seconds" << std::endl;

    // // 划分
    // // 哈希划分message分片，传入参数：哈希桶数目
    // // divideMessage(100, comment_hasCreator_person_r, post_hasCreator_person_r,comment_hasCreator_person, post_hasCreator_person,comment_0_0,post_0_0,comment_replyOf_comment_r,comment_replyOf_post_r,comment_replyOf_comment,comment_replyOf_post,comment_hasTag_tag,post_hasTag_tag,comment_isLocatedIn_place,post_isLocatedIn_place);
    
    // // divideMessageHashAll(100, comment_hasCreator_person_r, post_hasCreator_person_r,comment_hasCreator_person, post_hasCreator_person,forum_containerOf_post_r,comment_0_0,post_0_0,forum_0_0,comment_replyOf_comment_r,comment_replyOf_post_r,comment_replyOf_comment,comment_replyOf_post,comment_hasTag_tag,post_hasTag_tag,comment_isLocatedIn_place,post_isLocatedIn_place);
    
    // // sf0.1
    // divideMessageCommunity(10, 286744, person_0_0, person_likes_comment_r, person_likes_post_r, comment_hasCreator_person_r, post_hasCreator_person_r,comment_hasCreator_person, post_hasCreator_person,forum_containerOf_post_r,comment_0_0,post_0_0,comment_replyOf_comment_r,comment_replyOf_post_r,comment_replyOf_comment,comment_replyOf_post,comment_hasTag_tag,post_hasTag_tag,comment_isLocatedIn_place,post_isLocatedIn_place);
    // // divideMessageCommunityThreads(10, 286744, person_0_0, person_likes_comment_r, person_likes_post_r, comment_hasCreator_person_r, post_hasCreator_person_r,comment_hasCreator_person, post_hasCreator_person,forum_containerOf_post_r,comment_0_0,post_0_0,forum_0_0,comment_replyOf_comment_r,comment_replyOf_post_r,comment_replyOf_comment,comment_replyOf_post,comment_hasTag_tag,post_hasTag_tag,comment_isLocatedIn_place,post_isLocatedIn_place);
    // // sf3
    // // divideMessageCommunity(100, 9010236, person_0_0, person_likes_comment_r, person_likes_post_r, comment_hasCreator_person_r, post_hasCreator_person_r,comment_hasCreator_person, post_hasCreator_person,forum_containerOf_post_r,comment_0_0,post_0_0,forum_0_0,comment_replyOf_comment_r,comment_replyOf_post_r,comment_replyOf_comment,comment_replyOf_post,comment_hasTag_tag,post_hasTag_tag,comment_isLocatedIn_place,post_isLocatedIn_place);
    // // divideMessageCommunityThreads(10, 9010236, person_0_0, person_likes_comment_r, person_likes_post_r, comment_hasCreator_person_r, post_hasCreator_person_r,comment_hasCreator_person, post_hasCreator_person,forum_containerOf_post_r,comment_0_0,post_0_0,forum_0_0,comment_replyOf_comment_r,comment_replyOf_post_r,comment_replyOf_comment,comment_replyOf_post,comment_hasTag_tag,post_hasTag_tag,comment_isLocatedIn_place,post_isLocatedIn_place);
    // // sf30
    // // divideMessageCommunity(100, 87095182, person_0_0, person_likes_comment_r, person_likes_post_r, comment_hasCreator_person_r, post_hasCreator_person_r,comment_hasCreator_person, post_hasCreator_person,forum_containerOf_post_r,comment_0_0,post_0_0,forum_0_0,comment_replyOf_comment_r,comment_replyOf_post_r,comment_replyOf_comment,comment_replyOf_post,comment_hasTag_tag,post_hasTag_tag,comment_isLocatedIn_place,post_isLocatedIn_place);
    // // divideMessageCommunityThreads(10, 87095182, person_0_0, person_likes_comment_r, person_likes_post_r, comment_hasCreator_person_r, post_hasCreator_person_r,comment_hasCreator_person, post_hasCreator_person,forum_containerOf_post_r,comment_0_0,post_0_0,forum_0_0,comment_replyOf_comment_r,comment_replyOf_post_r,comment_replyOf_comment,comment_replyOf_post,comment_hasTag_tag,post_hasTag_tag,comment_isLocatedIn_place,post_isLocatedIn_place);
    // GetCurrentProcessMemoryUsage();
    // personAdd(comment_0_0,post_0_0);
    // forumAdd(post_0_0);
    // // forumAdd(forum_containerOf_post_r,post_hasCreator_person);

    // version2: 分步载入划分，降低内存占用
    loadMessageEdge(initialDir,person_0_0,comment_hasCreator_person_r,post_hasCreator_person_r,comment_hasCreator_person,post_hasCreator_person,post_hasTag_tag,forum_containerOf_post_r,comment_replyOf_comment_r,comment_replyOf_post_r,comment_replyOf_comment,comment_replyOf_post);
    loadMessageCreationDate(initialDir,comment_creationDate,post_creationDate);
    createLGraph(person_0_0,comment_hasCreator_person_r,post_hasCreator_person_r,comment_hasCreator_person,post_hasCreator_person,comment_replyOf_comment,comment_replyOf_post,g);
    std::cout<<"step1 start memory:"<<std::endl;
    GetCurrentProcessMemoryUsage();
    // divideCommunity(part, 286744, g, person_0_0, comment_hasCreator_person_r, post_hasCreator_person_r,comment_replyOf_comment,comment_replyOf_post,communitiesSet,commentSetArray,postSetArray);
    // divideCommunity(part, 9010236, g, person_0_0, comment_hasCreator_person_r, post_hasCreator_person_r,comment_replyOf_comment,comment_replyOf_post,communitiesSet,commentSetArray,postSetArray);
    divideCommunity(part, 87095182, g, person_0_0, comment_hasCreator_person_r, post_hasCreator_person_r,comment_replyOf_comment,comment_replyOf_post,communitiesSet,commentSetArray,postSetArray);
    divideMessageChain(communitiesSet,commentSetArray,postSetArray,comment_hasCreator_person,post_hasCreator_person,forum_containerOf_post_r,comment_replyOf_comment_r,comment_replyOf_post_r,comment_replyOf_comment,comment_replyOf_post,post_hasTag_tag,comment_creationDate,post_creationDate);
    std::cout<<"step1 end memory:"<<std::endl;
    GetCurrentProcessMemoryUsage();
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
    communitiesSet.clear();
    releaseNestedMap(communitiesSet);
    std::cout<<"step2 start memory:"<<std::endl;
    GetCurrentProcessMemoryUsage();
    loadMessageProperty(initialDir, comment_0_0, post_0_0,person_likes_comment_r,person_likes_post_r,post_isLocatedIn_place,comment_hasTag_tag,comment_isLocatedIn_place); 
    divideMessageProperty(commentSetArray,postSetArray,comment_0_0,post_0_0,comment_hasTag_tag,comment_isLocatedIn_place,post_isLocatedIn_place,person_likes_comment_r,person_likes_post_r);
    personAdd(comment_0_0,post_0_0);
    forumAdd(post_0_0);
    auto loadEnd = std::chrono::high_resolution_clock::now();

    // 获取结束时间点
    auto end = std::chrono::high_resolution_clock::now();

    auto outputDuration = std::chrono::duration_cast<std::chrono::seconds>(end - loadEnd);
    // std::cout << "write files uses: " << outputDuration.count() << " seconds" << std::endl;

    // 计算所消耗的时间
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << "Total time: " << duration.count() << " seconds" << std::endl;


    return 0;
}
