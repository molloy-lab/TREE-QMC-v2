#ifndef TREE_HPP
#define TREE_HPP

#include "utility.hpp"
#include "dict.hpp"
#include "taxa.hpp"

class Node {
    friend class Tree;
    friend class SpeciesTree;
    public:
        Node(index_t index);
        ~Node();
        void new_states(index_t size);
        void delete_states();
    private:
        Node *parent;
        std::vector<Node *> children;
        index_t index, size;
        weight_t /* **doublet, */ *singlet, s1, s2;
        // std::map<index_t, weight_t> doublet;
        std::vector<std::pair<index_t, weight_t>> *doublet;
        static weight_t get_doublet(weight_t *singlet, weight_t s1, weight_t s2, index_t x, index_t y);
        weight_t get_doublet(index_t a, index_t b);
        void add_doublet(index_t a, index_t b, weight_t c);
};

class Tree {
    public:
        Tree();
        Tree(const std::string &newick, Dict *dict);
        virtual ~Tree();
        std::string to_string();
        void resolve();
        index_t size();
        std::unordered_map<index_t, index_t> &get_indices();
        weight_t ***build_graph(Taxa &subset);
    protected:
        Node *root;
        std::unordered_map<index_t, Node*> index2node;
        Dict *dict;
        index_t pseudonym();
        std::string display_tree(Node *root);
    private:
        index_t pseudonyms;
        std::unordered_map<index_t, index_t> indices;
        void clear_states(Node *root);
        void build_states(Node *root, Taxa &subset);
        void depth(Node *root, index_t depth);
        weight_t get_doublet(Node *subtree, index_t x, index_t y, bool complement);
        void sa_doublet(Node *root, weight_t sum, index_t x);
        weight_t aa_doublet(Node *root, index_t x, index_t y);
        static bool cmp(const std::pair<index_t, weight_t> &a, const std::pair<index_t, weight_t> &b);
        void sort_doublet(Node *root);
        std::unordered_set<index_t> bad_edges(Node *root, Taxa &subset, weight_t ***graph);
        void good_edges(Node *root, Taxa &subset, weight_t ***graph);
        Node *build_tree(const std::string &newick);
        Node *build_subtree_from(Node *root);
        void resolve_tree(Node *root);
        void add_indices(Node *root, std::vector<index_t> &indices);
};

class SpeciesTree : public Tree {
    public:
        SpeciesTree(std::vector<Tree *> &input, Dict *dict, std::string mode);
        ~SpeciesTree();
    private:
        std::string mode;
        Node *construct_stree(std::vector<Tree *> &input, Taxa &subset);
        Node *reroot(Node *root, std::unordered_set<index_t> &visited);
        Node *reroot_stree(Node *root, index_t pseudo);
        Node *pseudo2node(Node *root, index_t pseudo);
};

#endif
