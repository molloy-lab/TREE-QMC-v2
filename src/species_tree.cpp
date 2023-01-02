#include "tree.hpp"
#include "graph.hpp"

SpeciesTree::SpeciesTree(std::vector<Tree *> &input, Dict *dict, std::string mode) {
    this->dict = dict;
    this->artifinyms = dict->max_size();
    this->mode = mode;
    Taxa subset(dict, mode);
    switch (mode[1]) {
        case '0': {
            root = construct_stree(input, subset);
            break;
        }
        case '1': {
            std::unordered_map<quartet_t, weight_t> quartets;
            for (Tree * tree: input) {
                tree->get_quartets(&quartets);
            }
            // std::cout << to_string(quartets);
            // std::cout << quartets.size() << std::endl;
            root = construct_stree(quartets, subset);
            break;
        }
        case '2': {
            break;
        }
        default: {
            break;
        }
    }
    // std::cout << artifinyms << std::endl;
}

SpeciesTree::~SpeciesTree() {
    
}

index_t SpeciesTree::artifinym() {
    return -- artifinyms;
}

Node *SpeciesTree::construct_stree(std::vector<Tree *> &input, Taxa &subset) {
    index_t size = subset.size();
    Node *root;
    if (size < 4) {
        if (size == 1) {
            root = new Node(subset.root_at(0));
        }
        else if (size == 2) {
            root = new Node(pseudonym());
            root->children.push_back(new Node(subset.root_at(0)));
            root->children.push_back(new Node(subset.root_at(1)));
            root->children[0]->parent = root->children[1]->parent = root;
        }
        else {
            root = new Node(pseudonym());
            root->children.push_back(new Node(pseudonym()));
            Node *left = root->children[0];
            left->children.push_back(new Node(subset.root_at(0)));
            left->children.push_back(new Node(subset.root_at(1)));
            left->children[0]->parent = left->children[1]->parent = left;
            root->children.push_back(new Node(subset.root_at(2)));
            root->children[0]->parent = root->children[1]->parent = root;
        }
    }
    else {
        Graph *g = new Graph(input, subset);
        std::vector<index_t> A, B;
        weight_t max = g->get_cut(&A, &B);
        if (max < 0) {
            root = new Node(pseudonym());
            for (index_t i = 0; i < subset.size(); i ++) {
                Node *child = new Node(subset.root_at(i));
                root->children.push_back(child);
                child->parent = root;
            }
        }
        else {
            Taxa subsetA(subset), subsetB(subset);
            index_t artificial = artifinym();
            subsetA.struct_update(A, artificial);
            subsetB.struct_update(B, artificial);
            root = new Node(pseudonym());
            root->children.push_back(reroot_stree(construct_stree(input, subsetA), artificial));
            root->children.push_back(reroot_stree(construct_stree(input, subsetB), artificial));
            root->children[0]->parent = root->children[1]->parent = root;
        }
        delete g;
    }
    // std::cout << display_tree(root) << std::endl;
    return root;
}

Node *SpeciesTree::construct_stree(std::unordered_map<quartet_t, weight_t> &quartets, Taxa &subset) {
    index_t size = subset.size();
    Node *root;
    if (size < 4) {
        if (size == 1) {
            root = new Node(subset.root_at(0));
        }
        else if (size == 2) {
            root = new Node(pseudonym());
            root->children.push_back(new Node(subset.root_at(0)));
            root->children.push_back(new Node(subset.root_at(1)));
            root->children[0]->parent = root->children[1]->parent = root;
        }
        else {
            root = new Node(pseudonym());
            root->children.push_back(new Node(pseudonym()));
            Node *left = root->children[0];
            left->children.push_back(new Node(subset.root_at(0)));
            left->children.push_back(new Node(subset.root_at(1)));
            left->children[0]->parent = left->children[1]->parent = left;
            root->children.push_back(new Node(subset.root_at(2)));
            root->children[0]->parent = root->children[1]->parent = root;
        }
    }
    else {
        // std::cout << to_string(quartets) << std::endl;
        Graph *g = new Graph(quartets, subset);
        std::vector<index_t> A, B;
        g->get_cut(&A, &B);
        std::unordered_set<index_t> setA(A.begin(), A.end()), setB(B.begin(), B.end());
        Taxa subsetA(subset), subsetB(subset);
        index_t artificial = artifinym();
        subsetA.struct_update(B, artificial);
        subsetB.struct_update(A, artificial);
        root = new Node(pseudonym());
        std::unordered_map<quartet_t, weight_t> quartetsA, quartetsB;
        for (auto elem : quartets) {
            index_t *indices = split(elem.first);
            index_t count = 0;
            for (index_t i = 0; i < 4; i ++) {
                if (setA.find(indices[i]) != setA.end()) 
                    count ++;
            }
            switch (count) {
                case 0: {
                    if (quartetsB.find(elem.first) == quartetsB.end()) 
                        quartetsB[elem.first] = 0;
                    quartetsB[elem.first] += elem.second;
                    break;
                }
                case 1: {
                    for (index_t i = 0; i < 4; i ++) {
                        if (setA.find(indices[i]) != setA.end()) {
                            indices[i] = artificial;
                        }
                    }
                    quartet_t temp = join(indices);
                    if (quartetsB.find(temp) == quartetsB.end()) 
                        quartetsB[temp] = 0;
                    quartetsB[temp] += elem.second / A.size();
                    break;
                }
                case 2: {
                    break;
                }
                case 3: {
                    for (index_t i = 0; i < 4; i ++) {
                        if (setB.find(indices[i]) != setB.end()) {
                            indices[i] = artificial;
                        }
                    }
                    quartet_t temp = join(indices);
                    if (quartetsA.find(temp) == quartetsA.end()) 
                        quartetsA[temp] = 0;
                    quartetsA[temp] += elem.second / B.size();
                    break;
                }
                case 4: {
                    if (quartetsA.find(elem.first) == quartetsA.end()) 
                        quartetsA[elem.first] = 0;
                    quartetsA[elem.first] += elem.second;
                    break;
                }
            }
            delete [] indices;
        }
        root->children.push_back(reroot_stree(construct_stree(quartetsA, subsetA), artificial));
        root->children.push_back(reroot_stree(construct_stree(quartetsB, subsetB), artificial));
        root->children[0]->parent = root->children[1]->parent = root;
        delete g;
    }
    std::cout << display_tree(root) << std::endl;
    return root;
}

Node *SpeciesTree::reroot(Node *root, std::unordered_set<index_t> &visited) {
    std::vector<Node *> child;
    if (root->parent != NULL && visited.find(root->parent->index) == visited.end()) {
        visited.insert(root->parent->index);
        child.push_back(reroot(root->parent, visited));
    }
    for (Node *ch : root->children) {
        if (ch != NULL && visited.find(ch->index) == visited.end()) {
            visited.insert(ch->index);
            child.push_back(reroot(ch, visited));
        }
    }
    Node *new_root;
    if (child.size() >= 2) {
        new_root = new Node(pseudonym());
        visited.insert(new_root->index);
        for (Node *ch : child) {
            new_root->children.push_back(ch);
            ch->parent = new_root;
        }
        return new_root;
    }
    else if (child.size() == 1) {
        new_root = child[0];
    }
    else {
        new_root = new Node(root->index);
    }
    // std::cout << '>' << display_tree(new_root) << std::endl;
    return new_root;
}

Node *SpeciesTree::reroot_stree(Node *root, index_t artificial) {
    Node *new_root = artificial2node(root, artificial);
    std::unordered_set<index_t> visited;
    visited.insert(new_root->index);
    Node *new_tree = reroot(new_root, visited);
    delete root;
    return new_tree;
}

Node *SpeciesTree::artificial2node(Node *root, index_t artificial) {
    if (root->children.size() == 0) {
        if (root->index == artificial) 
            return root;
        return NULL;
    }
    else {
        for (Node *child : root->children) {
            Node *temp = artificial2node(child, artificial);
            if (temp != NULL) return temp;
        }
        return NULL;
    }
}
