#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <cstring>
#include <vector>


template<typename ValueType>
class Trie
{
public:
    Trie();
    ~Trie();
    void reset();
    void insert(const std::string& key, const ValueType& value);
    std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const;

//    void dump();                    // remember to comment out
    
      // C++11 syntax for preventing copying and assignment
    Trie(const Trie&) = delete;
    Trie& operator=(const Trie&) = delete;
private:
    struct Node{
        char label;
        std::vector<Node*> children;
        std::vector<ValueType> values;
    };
    
    Node* m_root;
    void deleteTrie(Node* n);
    bool insertHelper(const char key[], const ValueType& value, Node* curr);
    std::vector<ValueType> findHelper(const char key[], bool exactMatchOnly, std::vector<ValueType>& matches, Node* curr) const;
    
//    void toilet(Node* n);
};

template<typename ValueType>
Trie<ValueType>::Trie(){
    m_root = new Node();
}

template<typename ValueType>
Trie<ValueType>::~Trie(){
    deleteTrie(m_root);
}

template<typename ValueType>
void Trie<ValueType>::reset(){
    deleteTrie(m_root);
    m_root = new Node();
}

template<typename ValueType>
void Trie<ValueType>::deleteTrie(Node* n){     // deletes all nodes in Trie
    if(n == nullptr)
        return;
    
    for(int i = 0; i < n->children.size(); i++){
        deleteTrie(n->children[i]);
    }
    delete n; 
}

template<typename ValueType>
void Trie<ValueType>::insert(const std::string& key, const ValueType& value){
    if(key.length() == 0)   // check that key is valid (is not empty)
        return;
    
    insertHelper(key.c_str(), value, m_root);
}

template<typename ValueType>
bool Trie<ValueType>::insertHelper(const char key[], const ValueType& value, Node* curr){
    if(key[0] == '\0')      // returns true if the key is empty, returns false otherwise
        return true;

    for(int i = 0; i < curr->children.size(); i++){   // check if the next label exists already
        if(curr->children[i]->label == key[0]){
            // if the next key char is null, then add the value to the current node
            if(insertHelper(key+1, value, curr->children[i])){
                curr->children[i]->values.push_back(value);
            }
            return false;
        }
    }
    
    Node* n = new Node();
    n->label = key[0];
    curr->children.push_back(n);
    if(insertHelper(key+1, value, n)){
        n->values.push_back(value);
    }
    return false;
}

template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::find(const std::string& key, bool exactMatchOnly) const{
    std::vector<ValueType> matches;
    Node* curr = m_root;
    
    // check that the first char is a match, regardless of exact matches
    for(int i = 0; i < curr->children.size(); i++){
        if(curr->children[i] != nullptr && curr->children[i]->label == key[0]){
            return findHelper(key.c_str(), exactMatchOnly, matches, curr->children[i]);
        }
    }
    // if there are no matching first chars, return an empty vector
    return matches;
}

template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::findHelper(const char key[], bool exactMatchOnly, std::vector<ValueType>& matches, Node* curr) const{
    
    // base case: key is empty, return immediately
    if(key[0] == '\0')
        return matches;
    
    // base case 2: the node is a nullptr, return immediately
    if(curr == nullptr)
        return matches;
    
    // the label of the current node equals the next key char
    if(curr->label == key[0]){
        // if the key char is the last one, add all values in the node to the vector
        if(key[1] == '\0'){
            for(int i = 0; i < curr->values.size(); i++){
                matches.push_back(curr->values[i]);
            }
        }else{
            // if the key char is not the last one, call findHelper() on all of the curr's children
            for(int i = 0; i < curr->children.size(); i++){
                findHelper(key+1, exactMatchOnly, matches, curr->children[i]);
            }
        }
    // if we are not looking for exact matches, call findHelper() on curr's children
    }else if(!exactMatchOnly){
        if(key[1] == '\0'){
            for(int i = 0; i < curr->values.size(); i++){
                matches.push_back(curr->values[i]);
            }
        }else{
            for(int i = 0; i < curr->children.size(); i++){
                findHelper(key+1, true, matches, curr->children[i]);
            }
        }
    }
    return matches;
}



//////////////////////////////////
/*
template<typename ValueType>
void Trie<ValueType>::dump(){
    toilet(m_root);
}

template<typename ValueType>
void Trie<ValueType>::toilet(Node* n){
    if(n == nullptr )
        return;
    
    for(int i = 0; i < n->values.size(); i++){
        std::cout << n->label << "\t\t" << n->values[i] << std::endl;
    }
    
    for(int i = 0; i < n->children.size(); i++){
        toilet(n->children[i]);
    }
}
*/
#endif // TRIE_INCLUDED
