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

#include <sys/stat.h>


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

    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> creator1;
    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> creator2;
    
    
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
    file2.open(directory1+"/comment_hasCreator_person_0_0.csv");
    file2<<"Comment.id|Person.id"<<std::endl;

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
            file2<<line<<std::endl;
        }
    }
    file.close();
    file2.close();

    // post hasCreator person
    file.open(directory+"/post_hasCreator_person_0_0.csv");
    file2.open(directory1+"/post_hasCreator_person_0_0.csv");
    file2<<"Post.id|Person.id"<<std::endl;

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
    for(auto i=commentPath.begin();i!=commentPath.end();i++) {
        std::string from = i->first;
    	std::vector<std::string> tmp = i->second;
        if(tmp.size()==1) {
            // comment没有comment邻居
            // ！！！这里也要存post
            comment.insert(tmp[0]);
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
                // if(visited1[tmp[2]].find(j)==visited1[tmp[2]].end()) {
                //     visited1[tmp[2]][j] = 1;
                //     file3_2<<tmp[2]<<'|'<<j<<std::endl;
                // }
                // file3_2<<tmp[1]<<'|'<<j<<std::endl;
            }
        }
    }
	
	// comment_hasCreator_person
    file.open(directory+"/comment_hasCreator_person_0_0.csv");
    file2.open(directory1+"/comment_hasCreator_person_0_0.csv");
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
    file2.open(directory1+"/post_hasCreator_person_0_0.csv");
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

    
	// forum_hasTag_tag
    file.open(directory+"/forum_hasTag_tag_0_0.csv");
    file2.open(directory1+"/forum_hasTag_tag_0_0.csv");
    
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
        
        if(forum.find(items[0])!=forum.end()){
            file2<<line<<std::endl;
            // forum.insert(items[0]); // 这里加不加都无所谓，提前已经完成了person_forum
            // newForum.insert(items[0]);
        }
    }
    file.close();
    file2.close(); 

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

int main() {
    // 获取开始时间点
    auto start = std::chrono::high_resolution_clock::now();

    std::ifstream file("./input/social_network-csv_composite-longdateformatter-sf0.1/dynamic/person_0_0.csv");
    std::string line;
    std::map<std::string, Person> people;
    std::map<std::string, std::map<std::string,std::string>> knows;

    // Skip header line
    // std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> items;

        while (std::getline(ss, item, '|')) {
            items.push_back(item);
        }

        if (items.size() != 10) {
            std::cerr << "Invalid line: " << line << std::endl;
            continue;
        }

        Person person;
        person.id = items[0];
        person.line = line;

        people[items[0]]=person;
    }

    file.close();
    std::cout<<"all Person num:"<<people.size()<<std::endl;

    // special
    // 

    // Print information for each person
    // for (const auto& p : people) {
    //     p.printInfo();
    // }


    // Read data from person_knows_person_0_0.csv
    std::ifstream file2("./input/social_network-csv_composite-longdateformatter-sf0.1/dynamic/person_knows_person_0_0.csv");
    // std::getline(file2, line); // Skip header line

    while (std::getline(file2, line)) {
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

        std::string source = items[0];
        std::string target = items[1];
        std::string content = items[2];

        knows[source][target]=content;
        knows[target][source]=content;
    }

    file2.close();

    // 计算所有person的一跳邻居
    std::map<std::string, int> neighbor_counts;
    calculateNeighborCounts(people, knows, neighbor_counts);

    
    // 计算所有person的两跳邻居
    std::map<std::string, std::set<std::string>> two_hop_neighbors;
    std::map<std::string, int> two_hop_neighbor_counts = calculateTwoHopNeighborCounts(people, knows, two_hop_neighbors);

    // 计算所有person的一跳邻居数目平均数
    double average_neighbor_count = 0.0;
    for (const auto& pair : neighbor_counts) {
        average_neighbor_count += pair.second;
    }
    average_neighbor_count /= neighbor_counts.size();
    std::cout << "Average number of one-hop neighbors: " << average_neighbor_count << std::endl;

    // 计算所有person的两跳邻居数目平均数
    double average_two_hop_neighbor_count = 0.0;
    for (const auto& pair : two_hop_neighbor_counts) {
        average_two_hop_neighbor_count += pair.second;
    }
    average_two_hop_neighbor_count /= neighbor_counts.size();
    std::cout << "Average number of two-hop neighbors: " << average_two_hop_neighbor_count << std::endl;
    
    // 划分方法1：每个person做中心person
    int part = 0;
    for(auto it=people.begin();it!=people.end();it++){
        std::cout<<part<<std::endl;
        // if(part>2) break;
        
        std::string givenPersonID = it->first; // Example ID
        std::set<std::string> twoHopResults;

        // getKnows(givenPersonID, knows, twoHopResults);
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

        getDivide(twoHopResults, "./input/social_network-csv_composite-longdateformatter-sf0.1/dynamic/", "./output/"+givenPersonID+"/", comment, post, knows);
        
        std::cout<<"part persons num:"<<twoHopResults.size()<<std::endl<<std::endl;
        part++;
    }

    // 划分方法2：贪心算法合并person
    // int part = 0;
    // std::set<std::string> peopleSet;
    // for(auto it=people.begin();it!=people.end();it++){
    //     peopleSet.insert(it->first);
    // }
    // std::ofstream tmpfile("./output/partTwoHopPersons.txt");
    
    // while(peopleSet.size()!=0) {
    //     std::string seed = getBiggestOneHop(peopleSet,neighbor_counts);
    //     if(seed=="") continue;
    //     // std::cout<<seed<<std::endl;
    //     std::set<std::string> twoHopResults;
    //     twoHopResults.insert(seed);
    //     twoHopResults.insert(two_hop_neighbors[seed].begin(), two_hop_neighbors[seed].end());
    //     // std::cout<<"twohopsize:"<<twoHopResults.size()<<std::endl;
    //     std::set<std::string> center;
    //     center.insert(seed);
    //     std::set<std::string> cNeighbor = getKnowsNeighbor1(center,knows,peopleSet);
    //     // int threshhold = calculateSetTwoHopNeighborNum(center, two_hop_neighbor_counts);
    //     int threshhold = twoHopResults.size();
    //     for (const auto& person : center) {
    //         cNeighbor.erase(person);
    //     }
    //     // while(threshhold<100000&&cNeighbor.size()!=0) {
    //     while(cNeighbor.size()!=0) {
    //         // std::cout<<"cNeighborsize1:"<<cNeighbor.size()<<std::endl;
    //         std::string smallestOneHop = getSmallestOneHop(cNeighbor,neighbor_counts);
    //         if(smallestOneHop=="") break;
    //         threshhold = calculateSetTwoHopNeighborNumNew(threshhold, smallestOneHop, twoHopResults, two_hop_neighbors); 
    //         if(threshhold>1700) break;
    //         // std::cout<<"smallest onehop:"<<smallestOneHop<<std::endl;
    //         center.insert(smallestOneHop);
    //         twoHopResults.insert(smallestOneHop);
    //         twoHopResults.insert(two_hop_neighbors[smallestOneHop].begin(), two_hop_neighbors[smallestOneHop].end());
            
    //         // threshhold = calculateSetTwoHopNeighborNum(center, two_hop_neighbor_counts); 
            
    //         // std::cout<<"threshhold:"<<threshhold<<std::endl;
    //         // 获取center的邻居
    //         cNeighbor = getKnowsNeighbor1(center,knows,peopleSet);
    //         // std::cout<<"cNeighborsize:"<<cNeighbor.size()<<std::endl;
    //         for (const auto& person : center) {
    //             cNeighbor.erase(person);
    //         }
            
    //     }
    //     std::cout<<"center persons num:"<<center.size()<<std::endl;
    //     std::cout<<"threshhold:"<<threshhold<<std::endl;
    //     std::cout<<"twohop persons num:"<<twoHopResults.size()<<std::endl;
    //     std::string outPath = "./output/";
    //     outPath+=std::to_string(part);
    //     for(auto it=center.begin();it!=center.end();it++){
    //         tmpfile<<part<<' '<<*it<<std::endl;
    //     }

    //     std::set<std::string> comment;
    //     std::set<std::string> post;
    //     getDivide(twoHopResults, "./input/social_network-csv_composite-longdateformatter-sf0.1/dynamic/", outPath, comment, post, knows);
        
    //     std::cout<<"part persons num:"<<twoHopResults.size()<<std::endl<<std::endl;
        
    //     part++;
    //     std::cout<<"part:"<<part<<std::endl;

    //     // peopleSet.erase(center.begin(),center.end());
    //     for (const auto& person : center) {
    //         peopleSet.erase(person);
    //     }
    // }

    // 划分方法3：3.1贪心算法合并person，3.2合并相似度高的两跳person集
    // std::set<std::set<std::string>> twoHopSets;
    // std::set<std::string> peopleSet;
    // for(auto it=people.begin();it!=people.end();it++){
    //     peopleSet.insert(it->first);
    // }
    // int testCenter=0;
    // int testtwohop=0;
    // while(peopleSet.size()!=0) {
    //     // std::string seed = getSmallestOneHop(peopleSet,neighbor_counts);
    //     std::string seed = getBiggestOneHop(peopleSet,neighbor_counts);
    //     if(seed=="") continue;
    //     // std::cout<<seed<<std::endl;
    //     std::set<std::string> twoHopResults;
    //     twoHopResults.insert(seed);
    //     twoHopResults.insert(two_hop_neighbors[seed].begin(), two_hop_neighbors[seed].end());
    //     // std::cout<<"twohopsize:"<<twoHopResults.size()<<std::endl;
    //     std::set<std::string> center;
    //     center.insert(seed);
    //     std::set<std::string> cNeighbor = getKnowsNeighbor1(center,knows,peopleSet);
    //     // int threshhold = calculateSetTwoHopNeighborNum(center, two_hop_neighbor_counts);
    //     int threshhold = twoHopResults.size();
    //     for (const auto& person : center) {
    //         cNeighbor.erase(person);
    //     }
    //     while(cNeighbor.size()!=0) {
    //         // std::cout<<"cNeighborsize1:"<<cNeighbor.size()<<std::endl;
    //         std::string smallestOneHop = getSmallestOneHop(cNeighbor,neighbor_counts);
    //         if (smallestOneHop == "") {
    //             break;
    //         }

    //         threshhold = calculateSetTwoHopNeighborNumNew(threshhold, smallestOneHop, twoHopResults, two_hop_neighbors);
    //         if(threshhold>110) break;

    //         // std::cout<<"smallest onehop:"<<smallestOneHop<<std::endl;
    //         center.insert(smallestOneHop);
    //         twoHopResults.insert(smallestOneHop);
    //         twoHopResults.insert(two_hop_neighbors[smallestOneHop].begin(), two_hop_neighbors[smallestOneHop].end());
            
    //         // threshhold = calculateSetTwoHopNeighborNum(center, two_hop_neighbor_counts); 
    //         // std::cout<<"threshhold:"<<threshhold<<std::endl;
    //         // 获取center的邻居
    //         cNeighbor = getKnowsNeighbor1(center,knows,peopleSet);
    //         // std::cout<<"cNeighborsize:"<<cNeighbor.size()<<std::endl;
    //         for (const auto& person : center) {
    //             cNeighbor.erase(person);
    //         }
            
    //     }
        
    //     // 插入到twoHopSets中待合并
    //     twoHopSets.insert(twoHopResults);
        

    //     for (const auto& person : center) {
    //         peopleSet.erase(person);
    //     }

    //     testCenter+=center.size();
    //     testtwohop+=twoHopResults.size();
    // }
    // std::cout<<"all center person num:"<<testCenter<<std::endl;
    // std::cout<<"all twohop person num:"<<testtwohop<<std::endl;
    
    // // 合并相似度高的两跳person集
    // double merge_threshold = 0.1;  // 相似度阈值
    // std::set<std::set<std::string>> mergedSets = merge_high_similarity_sets(twoHopSets, merge_threshold);
    // std::cout << "Number of merged sets: " << mergedSets.size() << std::endl;
    // std::cout << "Number of twoHopSets: " << twoHopSets.size() << std::endl;
    
    // /*twoHopSets.insert(peopleSet);
    // std::set<std::set<std::string>> mergedSets = twoHopSets;*/
    
    // int part = 0;
    // // 通过合并后的种子集合生成分片
    // std::ofstream tmpfile("./output/partTwoHopPersons.txt");
    // for(auto i: mergedSets) {
    //     std::string outPath = "./output/";
    //     outPath+=std::to_string(part);
    //     std::cout<<"part"<<part<<std::endl;
    //     std::cout<<"twohop persons num:"<<i.size()<<std::endl;
    //     for(auto it=i.begin();it!=i.end();it++){
    //         // outPath += *it;
    //         // outPath += "_";
    //         tmpfile<<part<<' '<<*it<<std::endl;
    //     }
    //     std::set<std::string> comment;
    //     std::set<std::string> post;
    //     getDivide(i, "./input/social_network-csv_composite-longdateformatter-sf0.1/dynamic/", outPath, comment, post, knows);
    //     std::cout<<"part persons num:"<<i.size()<<std::endl<<std::endl;
    //     part++;
    // }

    // 获取结束时间点
    auto end = std::chrono::high_resolution_clock::now();

    // 计算所消耗的时间
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << "Time taken by function: " << duration.count() << " seconds" << std::endl;


    return 0;
}
