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
        index_t index, size, depth;
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
        size_t resolve();
        index_t size();
        std::unordered_map<index_t, index_t> &get_indices();
        weight_t ***build_graph(Taxa &subset);
        void get_quartets(std::unordered_map<quartet_t, weight_t> *quartets);
        std::string to_string(std::unordered_map<quartet_t, weight_t> &quartets);
    protected:
        Node *root;
        std::unordered_map<index_t, Node*> index2node;
        Dict *dict;
        index_t pseudonym();
        std::string display_tree(Node *root);
        std::string display_tree_index(Node *root);
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
        size_t resolve_tree(Node *root);
        void add_indices(Node *root, std::vector<index_t> &indices);
        void get_leaves(Node *root, std::vector<Node *> *leaves);
        void get_depth(Node *root, index_t depth);
};

class SpeciesTree : public Tree {
    public:
        SpeciesTree(std::vector<Tree *> &input, Dict *dict, std::string mode);
        ~SpeciesTree();
    private:
        index_t artifinyms;
        std::string mode;
        index_t artifinym();
        Node *construct_stree(std::vector<Tree *> &input, Taxa &subset, index_t parent_pid, index_t depth);
        Node *construct_stree(std::unordered_map<quartet_t, weight_t> &input, Taxa &subset, index_t parent_pid, index_t depth);
        Node *reroot(Node *root, std::unordered_set<index_t> &visited);
        Node *reroot_stree(Node *root, index_t artificial);
        Node *artificial2node(Node *root, index_t artificial);
};

extern std::ofstream subproblem_csv;
extern std::string verbose;
extern unsigned long long count[8];

#endif
