#include "tree.hpp"
#include "graph.hpp"


SpeciesTree::SpeciesTree(std::vector<Tree *> &input, Dict *dict, std::string mode) {
    this->dict = dict;
    this->artifinyms = dict->size() * 2 - 3;
    this->mode = mode;
    Taxa subset(dict, mode[0]);
    root = construct_stree(input, subset);
    std::cout << artifinyms << std::endl;
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
        g->get_cut(&A, &B);
        Taxa subsetA(subset), subsetB(subset);
        index_t artificial = artifinym();
        subsetA.struct_update(A, artificial);
        subsetB.struct_update(B, artificial);
        root = new Node(pseudonym());
        root->children.push_back(reroot_stree(construct_stree(input, subsetA), artificial));
        root->children.push_back(reroot_stree(construct_stree(input, subsetB), artificial));
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
    if (child.size() >= 2) {
        Node *new_root = new Node(pseudonym());
        visited.insert(new_root->index);
        for (Node *ch : child) {
            new_root->children.push_back(ch);
            ch->parent = new_root;
        }
        return new_root;
    }
    else if (child.size() == 1) {
        return child[0];
    }
    else {
        Node *new_root = new Node(root->index);
        return new_root;
    }
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
