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
#include <filesystem>

#include <dirent.h> // 用于目录操作
#include <sys/stat.h> // 用于检查文件状态


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

void loadExknows1(std::string filename,
                std::map<std::string, std::map<std::string,std::string>>& knows,
                std::unordered_map<std::string,std::unordered_set<std::string>>& person_hasInterest_tag,
                std::unordered_map<std::string,std::unordered_set<std::string>>& person_isLocatedIn_place,
                std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& person_studyAt_organisation,
                std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& person_workAt_organisation,
                std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& forum_hasMember_person_r,
                std::unordered_map<std::string,std::unordered_map<std::string,std::string>>& forum_hasMember_person,//new
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
                std::unordered_map<std::string,std::unordered_set<std::string>>& forum_containerOf_post,//new
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
        forum_hasMember_person[items[0]][items[1]] = items[2];
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
        forum_containerOf_post[items[0]].insert(items[1]);
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

void checkEdgeValid(std::string inputdir,std::map<std::string, Person> people,std::map<std::string, std::map<std::string, std::string>>& knows) {
    // 划分方法3.1 一次load：3.1贪心算法合并person，3.2合并相似度高的两跳person集
    std::unordered_map<std::string,std::unordered_set<std::string>> person_hasInterest_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> person_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_studyAt_organisation;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_workAt_organisation;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> forum_hasMember_person_r;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> forum_hasMember_person;//new
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
    std::unordered_map<std::string,std::unordered_set<std::string>> forum_containerOf_post;//new
    std::unordered_map<std::string, std::unordered_set<std::string>> reply1;
    std::unordered_map<std::string, std::unordered_set<std::string>> reply2;
    std::unordered_map<std::string, std::string> person_0_0;
    std::unordered_map<std::string, std::string> comment_0_0;
    std::unordered_map<std::string, std::string> post_0_0;
    std::unordered_map<std::string, std::string> forum_0_0;

    loadExknows1(inputdir,knows, person_hasInterest_tag,person_isLocatedIn_place,person_studyAt_organisation,person_workAt_organisation,forum_hasMember_person_r,forum_hasMember_person,forum_hasModerator_person_r,forum_hasModerator_person,comment_hasCreator_person_r,comment_hasCreator_person,post_hasCreator_person_r,post_hasCreator_person,post_hasTag_tag,post_isLocatedIn_place,comment_hasTag_tag,comment_isLocatedIn_place,forum_containerOf_post_r,forum_containerOf_post,reply1,reply2,person_0_0,comment_0_0,post_0_0,forum_0_0);

    // for(auto i:person_hasInterest_tag) {
    //     if(people.find(i.first)==people.end()) {
    //         std::cout<<"person_hasInterest_tag error:"<<i.first<<std::endl;
    //     }
    // }
    // for(auto i:person_isLocatedIn_place) {
    //     if(people.find(i.first)==people.end()) {
    //         std::cout<<"person_isLocatedIn_place error:"<<i.first<<std::endl;
    //     }
    // }
    // for(auto i:person_studyAt_organisation) {
    //     if(people.find(i.first)==people.end()) {
    //         std::cout<<"person_studyAt_organisation error:"<<i.first<<std::endl;
    //     }
    // }
    // for(auto i:person_workAt_organisation) {
    //     if(people.find(i.first)==people.end()) {
    //         std::cout<<"person_workAt_organisation error:"<<i.first<<std::endl;
    //     }
    // }
    for(auto i:comment_hasCreator_person_r) {
        if(people.find(i.first)==people.end()) {
            std::cout<<"comment_hasCreator_person_r error:"<<i.first<<std::endl;
        }
    }
    // for(auto i:forum_hasMember_person_r) {
    //     if(people.find(i.first)==people.end()) {
    //         std::cout<<"forum_hasMember_person_r error:"<<i.first<<std::endl;
    //     }
    // }
    // for(auto i:forum_hasModerator_person_r) {
    //     if(people.find(i.first)==people.end()) {
    //         std::cout<<"forum_hasModerator_person_r error:"<<i.first<<std::endl;
    //     }
    // }
    // for(auto i:post_hasCreator_person_r) {
    //     if(people.find(i.first)==people.end()) {
    //         std::cout<<"post_hasCreator_person_r error:"<<i.first<<std::endl;
    //     }
    // }
    // for(auto i:forum_hasModerator_person) {
    //     if(forum_0_0.find(i.first)==forum_0_0.end()) {
    //         std::cout<<"forum_hasModerator_person error:"<<i.first<<std::endl;
    //     }
    // }
    // for(auto i:forum_containerOf_post) {
    //     if(forum_0_0.find(i.first)==forum_0_0.end()) {
    //         std::cout<<"forum_containerOf_post error:"<<i.first<<std::endl;
    //     }
    // }
    // for(auto i:forum_containerOf_post_r) {
    //     if(post_0_0.find(i.first)==post_0_0.end()) {
    //         std::cout<<"forum_containerOf_post_r error:"<<i.first<<std::endl;
    //     }
    // }
    for(auto i:comment_hasCreator_person) {
        if(comment_0_0.find(i.first)==comment_0_0.end()) {
            std::cout<<"comment_hasCreator_person error:"<<i.first<<std::endl;
        }
    }
    for(auto i:post_hasCreator_person) {
        if(post_0_0.find(i.first)==post_0_0.end()) {
            std::cout<<"post_hasCreator_person error:"<<i.first<<std::endl;
        }
    }
    // for(auto i:post_hasTag_tag) {
    //     if(post_0_0.find(i.first)==post_0_0.end()) {
    //         std::cout<<"post_hasTag_tag error:"<<i.first<<std::endl;
    //     }
    // }
    // for(auto i:post_isLocatedIn_place) {
    //     if(post_0_0.find(i.first)==post_0_0.end()) {
    //         std::cout<<"post_isLocatedIn_place error:"<<i.first<<std::endl;
    //     }
    // }
    for(auto i:comment_hasTag_tag) {
        if(comment_0_0.find(i.first)==comment_0_0.end()) {
            std::cout<<"comment_hasTag_tag error:"<<i.first<<std::endl;
        }
    }  
    for(auto i:comment_isLocatedIn_place) {
        if(comment_0_0.find(i.first)==comment_0_0.end()) {
            std::cout<<"comment_isLocatedIn_place error:"<<i.first<<std::endl;
        }
    }
    for(auto i:reply1) {
        if(comment_0_0.find(i.first)==comment_0_0.end()) {
            std::cout<<"reply1 error:"<<i.first<<std::endl;
        }
        for(auto j:i.second) {
            if(comment_0_0.find(j)==comment_0_0.end()) {
                std::cout<<"reply1 c1 error:"<<j<<std::endl;
            }
        }
    }
    for(auto i:reply2) {
        if(comment_0_0.find(i.first)==comment_0_0.end()) {
            std::cout<<"reply2 error:"<<i.first<<std::endl;
        }
        for(auto j:i.second) {
            if(post_0_0.find(j)==post_0_0.end()) {
                std::cout<<"reply2 post error:"<<j<<std::endl;
            }
        }
    }

}

void getnewdata(std::string inputdir,std::map<std::string, Person> people,std::map<std::string, std::map<std::string, std::string>>& knows) {
    // 划分方法3.1 一次load：3.1贪心算法合并person，3.2合并相似度高的两跳person集
    std::unordered_map<std::string,std::unordered_set<std::string>> person_hasInterest_tag;
    std::unordered_map<std::string,std::unordered_set<std::string>> person_isLocatedIn_place;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_studyAt_organisation;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> person_workAt_organisation;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> forum_hasMember_person_r;
    std::unordered_map<std::string,std::unordered_map<std::string,std::string>> forum_hasMember_person;//new
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
    std::unordered_map<std::string,std::unordered_set<std::string>> forum_containerOf_post;//new
    std::unordered_map<std::string, std::unordered_set<std::string>> reply1;
    std::unordered_map<std::string, std::unordered_set<std::string>> reply2;
    std::unordered_map<std::string, std::string> person_0_0;
    std::unordered_map<std::string, std::string> comment_0_0;
    std::unordered_map<std::string, std::string> post_0_0;
    std::unordered_map<std::string, std::string> forum_0_0;

    loadExknows1(inputdir,knows, person_hasInterest_tag,person_isLocatedIn_place,person_studyAt_organisation,person_workAt_organisation,forum_hasMember_person_r,forum_hasMember_person,forum_hasModerator_person_r,forum_hasModerator_person,comment_hasCreator_person_r,comment_hasCreator_person,post_hasCreator_person_r,post_hasCreator_person,post_hasTag_tag,post_isLocatedIn_place,comment_hasTag_tag,comment_isLocatedIn_place,forum_containerOf_post_r,forum_containerOf_post,reply1,reply2,person_0_0,comment_0_0,post_0_0,forum_0_0);
    std::set<std::string> commentSet;
    for(auto i:comment_hasCreator_person) {
        commentSet.insert(i.first);
    }
    for(auto i:comment_hasTag_tag) {
        commentSet.insert(i.first);
    }  
    for(auto i:comment_isLocatedIn_place) {
        commentSet.insert(i.first);
    }
    for(auto i:reply1) {
        commentSet.insert(i.first);
        for(auto j:i.second) {
            commentSet.insert(j);
        }
    }
    for(auto i:reply2) {
        commentSet.insert(i.first);
    }

    std::ifstream file;
    std::ofstream file2,file3,file4;
    std::string line;
    // comment
    file.open("/home/shared/data/social_network-csv_composite-longdateformatter-sf0.1/dynamic/comment_0_0.csv");
    file2.open("./comment_0_0.csv");
    file3.open("./comment_0_0_except.csv");
    file4.open("./comment_0_0_all.csv");
    
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }
        if(commentSet.find(items[0])!=commentSet.end()){
            file2<<line<<std::endl;
        } else {
            file3<<line<<std::endl;
        }
        file4<<line<<std::endl;
    }
    file.close();
    file2.close(); 
    file3.close(); 
    file4.close(); 
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

std::vector<std::string> readFileAndSort(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    std::vector<std::string> lines;

    while (std::getline(file, line)) {
        lines.push_back(line);
    }

    std::sort(lines.begin(), lines.end());
    return lines;
}

// 读取文件并排序
void readFileSortAndPrint(const std::string& filename,const std::string& filename1) {
    std::ifstream file(filename);
    std::ofstream file1(filename1);
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return;
    }

    std::string line;
    std::vector<std::string> lines;

    // 读取文件的每一行
    while (std::getline(file, line)) {
        lines.push_back(line);
    }

    // 对读取到的行进行排序
    std::sort(lines.begin(), lines.end());

    // 输出排序后的行
    for (const auto& sortedLine : lines) {
        file1 << sortedLine << std::endl;
    }
}

void compareFiles(const std::string& filename1, const std::string& filename2) {
    auto lines1 = readFileAndSort(filename1);
    auto lines2 = readFileAndSort(filename2);
    std::cout<<"lines1 size:"<<lines1.size()<<std::endl;
    std::cout<<"lines2 size:"<<lines2.size()<<std::endl;

    size_t maxSize = std::max(lines1.size(), lines2.size());
    for (size_t i = 0; i < maxSize; ++i) {
        if (i >= lines1.size() || i >= lines2.size() || lines1[i] != lines2[i]) {
            std::cout << "Line " << i + 1 << " is different." << std::endl;
            if (i < lines1.size()) {
                std::cout << "File 1: " << lines1[i] << std::endl;
            }
            if (i < lines2.size()) {
                std::cout << "File 2: " << lines2[i] << std::endl;
            }
        }
    }
}

// 计算文件的字符数
size_t countFileChars(const std::string& filename) {
    std::ifstream file(filename);
    size_t charCount = 0;
    char c;

    while (file.get(c)) {
        ++charCount;
    }

    return charCount;
}

// 比较两个文件的字符数
void compareFileChars(const std::string& filename1, const std::string& filename2) {
    size_t charCount1 = countFileChars(filename1);
    size_t charCount2 = countFileChars(filename2);

    std::cout << "File 1 (" << filename1 << ") has " << charCount1 << " characters." << std::endl;
    std::cout << "File 2 (" << filename2 << ") has " << charCount2 << " characters." << std::endl;

    if (charCount1 == charCount2) {
        std::cout << "Both files have the same number of characters." << std::endl;
    } else {
        std::cout << "The files have a different number of characters." << std::endl;
    }
}

std::vector<std::string> getSubdirectories(const std::string& directoryPath) {
    std::vector<std::string> subdirectories;
    DIR* dir = opendir(directoryPath.c_str()); // 打开目录
    struct dirent* entry;

    if (dir == nullptr) {
        std::cerr << "给定路径不存在或不是一个目录。" << std::endl;
        return subdirectories; // 返回空向量
    }

    while ((entry = readdir(dir)) != nullptr) { // 遍历目录项
        if (entry->d_type == DT_DIR) { // 检查是否为目录
            std::string name = entry->d_name;
            if (name == "." || name == "..") continue; // 跳过当前目录和上级目录
            subdirectories.push_back(directoryPath + name);
        }
    }
    closedir(dir); // 关闭目录
    return subdirectories;
}

int main() {
    // 获取开始时间点
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::string> files=getSubdirectories("/home/huzheyuan/dev/divide/output/");
    std::cout<<files.size()<<std::endl;
    for(auto i:files) {
        std::cout<<i+"----sort..."<<std::endl;
        readFileSortAndPrint(i+"/comment_0_0.csv",i+"/comment_0_0_sorted.csv");
        readFileSortAndPrint(i+"/forum_0_0.csv",i+"/forum_0_0_sorted.csv");
        readFileSortAndPrint(i+"/person_0_0.csv",i+"/person_0_0_sorted.csv");
        readFileSortAndPrint(i+"/post_0_0.csv",i+"/post_0_0_sorted.csv");
    
    }

    // readFileSortAndPrint("/home/huzheyuan/dev/gpstore_distributed/empty.csv","/home/huzheyuan/dev/gpstore_distributed/empty_sorted.csv");
    
    
    // compareFiles("/home/huzheyuan/dev/gpstore_distributed/comment_0_0_1.csv","/home/huzheyuan/dev/gpstore_distributed/comment_0_0_all.csv");
    // compareFileChars("/home/huzheyuan/dev/gpstore_distributed/comment_0_0_1.csv","/home/huzheyuan/dev/gpstore_distributed/comment_0_0_all.csv");

    // std::string initialDir = "/home/huzheyuan/dev/gpstore_distributed/input/social_network-csv_composite-longdateformatter-sf0.1/dynamic/";
    // initialDir = "/home/huzheyuan/dev/divide/output/10995116278436/";

    // // Read data from person_0_0.csv
    // std::ifstream file(initialDir+"person_0_0.csv");
    // std::string line;
    // std::map<std::string, Person> people;
    // std::map<std::string, std::map<std::string,std::string>> knows;

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

    //     Person person;
    //     person.id = items[0];
    //     person.line = line;

    //     people[items[0]]=person;
    // }

    // file.close();
    // std::cout<<"all Person num:"<<people.size()<<std::endl;


    // // Read data from person_knows_person_0_0.csv
    // std::ifstream file2(initialDir+"person_knows_person_0_0.csv");

    // while (std::getline(file2, line)) {
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

    //     std::string source = items[0];
    //     std::string target = items[1];
    //     std::string content = items[2];

    //     knows[source][target]=content;
    //     knows[target][source]=content;
    // }

    // file2.close();

    // // 计算所有person的一跳邻居
    // std::map<std::string, int> neighbor_counts;
    // calculateNeighborCounts(people, knows, neighbor_counts);

    
    // // 计算所有person的两跳邻居
    // std::map<std::string, std::set<std::string>> two_hop_neighbors;
    // std::map<std::string, int> two_hop_neighbor_counts = calculateTwoHopNeighborCounts(people, knows, two_hop_neighbors);

    // // 计算所有person的一跳邻居数目平均数
    // double average_neighbor_count = 0.0;
    // for (const auto& pair : neighbor_counts) {
    //     average_neighbor_count += pair.second;
    // }
    // average_neighbor_count /= neighbor_counts.size();
    // std::cout << "Average number of one-hop neighbors: " << average_neighbor_count << std::endl;

    // // 计算所有person的两跳邻居数目平均数
    // double average_two_hop_neighbor_count = 0.0;
    // for (const auto& pair : two_hop_neighbor_counts) {
    //     average_two_hop_neighbor_count += pair.second;
    // }
    // average_two_hop_neighbor_count /= neighbor_counts.size();
    // std::cout << "Average number of two-hop neighbors: " << average_two_hop_neighbor_count << std::endl;
    
    // // 划分方法1：每个person做中心person
    // // partition1("/home/huzheyuan/dev/gpstore_distributed/input/social_network-csv_composite-longdateformatter-sf0.1/dynamic/", people, knows, two_hop_neighbors);
    
    // // 划分方法1.1：一次load
    // // partition1_1("/home/huzheyuan/dev/gpstore_distributed/input/social_network-csv_composite-longdateformatter-sf0.1/dynamic/", people, knows, two_hop_neighbors);
    
    // // partition2_1(1200, "/home/huzheyuan/dev/gpstore_distributed/input/social_network-csv_composite-longdateformatter-sf0.1/dynamic/", people, knows, two_hop_neighbors, neighbor_counts);
    
    // // partition3_1(1100, 0.99, "/home/huzheyuan/dev/gpstore_distributed/input/social_network-csv_composite-longdateformatter-sf0.1/dynamic/", people, knows, two_hop_neighbors, neighbor_counts);

    // // checkEdgeValid(initialDir, people, knows);
    // getnewdata(initialDir, people, knows);

    // 获取结束时间点
    auto end = std::chrono::high_resolution_clock::now();

    // 计算所消耗的时间
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << "Time taken by function: " << duration.count() << " seconds" << std::endl;


    return 0;
}
