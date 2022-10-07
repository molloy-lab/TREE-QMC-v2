#include "tree.hpp"

extern long long total_count[2];

Tree::Tree() {
    pseudonyms = 0;
}

Tree::Tree(const std::string &newick, Dict *dict) {
    pseudonyms = 0;
    this->dict = dict;
    root = build_tree(newick);
}

Tree::~Tree() {
    delete root;
}

std::string Tree::to_string() {
    return display_tree(root) + ";";
}

size_t Tree::resolve() {
    return resolve_tree(root);
}

index_t Tree::size() {
    return index2node.size();
}

std::unordered_map<index_t, index_t> &Tree::get_indices() {
    return indices;
}

weight_t ***Tree::build_graph(Taxa &subset) {
    index_t s = subset.singleton_taxa(), m = subset.artificial_taxa();
    build_states(root, subset);
    for (index_t i = 0; i <= m; i ++) 
        sa_doublet(root, 0, i);
    for (index_t i = 1; i <= m; i ++) {
        for (index_t j = 0; j <= m; j ++) {
            if (i == j) continue;
            aa_doublet(root, i, j);
        }
    }
    weight_t ***graph = new weight_t**[2];
    graph[0] = Matrix::new_mat(subset.size());
    graph[1] = Matrix::new_mat(subset.size());
    bad_edges(root, subset, graph);
    good_edges(root, subset, graph);
    clear_states(root);
    return graph;
}

index_t Tree::pseudonym() {
    return - (++ pseudonyms);
}

Node *Tree::build_tree(const std::string &newick) {
    if (newick.length() == 0 || newick.at(0) != '(') {
        std::string label = newick.substr(0, newick.find(":"));
        //std::string length = newick.substr(newick.find(":") + 1, std::string::npos);
        Node *root = new Node(dict->label2index(label));
        index2node[root->index] = root;
        if (indices.find(root->index) == indices.end())
            indices[root->index] = 0;
        indices[root->index] ++;
        return root;
    }
    else {
        Node *root = new Node(pseudonym());
        int k = 1;
        for (int i = 0, j = 0; i < newick.length(); i ++) {
            if (newick.at(i) == '(') j ++;
            if (newick.at(i) == ')') j --;
            if (newick.at(i) == ',' && j == 1) {
                root->children.push_back(build_tree(newick.substr(k, i - k)));
                k = i + 1;
            }
        }
        int i = newick.length() - 1;
        while (newick.at(i) != ')') i --;
        /*
        std::string branch = newick.substr(i + 1, std::string::npos);
        if (branch.find(":") != std::string::npos) {
            std::string support = branch.substr(0, branch.find(":"));
            std::string length = branch.substr(branch.find(":") + 1, std::string::npos);
        }
        */
        root->children.push_back(build_tree(newick.substr(k, i - k)));
        for (Node *child : root->children) 
            child->parent = root;
        return root;
    }
}

std::string Tree::display_tree(Node *root) {
    if (root->children.size() == 0) 
        return dict->index2label(root->index);
    std::string s = "(";
    for (Node * node : root->children) 
        s += display_tree(node) + ",";
    s[s.size() - 1] = ')';
    return s;
}

std::string Tree::display_tree_index(Node *root) {
    if (root->children.size() == 0) 
        return std::to_string(root->index);
    std::string s = "(";
    for (Node * node : root->children) 
        s += display_tree_index(node) + ",";
    s[s.size() - 1] = ')';
    return s;
}

size_t Tree::resolve_tree(Node *root) {
    size_t total = 0;
    for (index_t i = 0; i < root->children.size(); i ++) 
        total += resolve_tree(root->children[i]);
    while (root->children.size() > 2) {
        index_t i = rand() % root->children.size(), j = i;
        while (j == i) j = rand() % root->children.size();
        Node *new_root = new Node(pseudonym());
        total ++;
        index2node[new_root->index] = new_root;
        new_root->children.push_back(root->children[i]);
        new_root->children.push_back(root->children[j]);
        root->children[i]->parent = root->children[j]->parent = new_root;
        root->children.erase(root->children.begin() + i);
        root->children.erase(root->children.begin() + (j > i ? j - 1 : j));
        root->children.push_back(new_root);
        new_root->parent = root;
    }
    return total;
}

void Tree::clear_states(Node *root) {
    root->delete_states();
    for (Node *child : root->children) 
        clear_states(child);
}

void Tree::build_states(Node *root, Taxa &subset) { 
    root->new_states(subset.artificial_taxa());
    if (root->children.size() == 0) {
        index_t index = subset.root_key(root->index);
        root->singlet[index] = subset.root_weight(root->index);
    }
    else {
        for (Node *child : root->children) {
            build_states(child, subset);
            for (index_t i = 0; i < root->size + 1; i ++) 
                root->singlet[i] += child->singlet[i];
        }
    }
    root->s1 = root->s2 = 0;
    for (int i = 1; i < root->size + 1; i ++) {
        root->s1 += root->singlet[i];
        root->s2 += root->singlet[i] * root->singlet[i];
    }
}

weight_t Tree::get_doublet(Node *subtree, index_t x, index_t y, bool complement) {
    if (complement) {
        weight_t *c = new weight_t[subtree->size + 1], s1 = 0, s2 = 0;
        for (int i = 0; i <= subtree->size; i ++) {
            c[i] = root->singlet[i] - subtree->singlet[i];
            if (i != 0) {
                s1 += c[i];
                s2 += c[i] * c[i];
            }
        }
        weight_t ret = Node::get_doublet(c, s1, s2, x, y);
        delete [] c;
        return ret;
    }
    return Node::get_doublet(subtree->singlet, subtree->s1, subtree->s2, x, y);
}

void Tree::sa_doublet(Node *root, weight_t sum, index_t x) {
    root->add_doublet(0, x, sum);
    // root->doublet[0][x] = sum;
    if (root->children.size() == 0) return ;
    sa_doublet(root->children[0], sum + get_doublet(root->children[1], x, 0, false), x);
    sa_doublet(root->children[1], sum + get_doublet(root->children[0], x, 0, false), x);
}

/*
weight_t Tree::aa_doublet(Node *root, index_t x, index_t y) {
    if (root->children.size() == 0) return 0;
    weight_t t0 = get_doublet(root->children[0], x, y, false), t1 = get_doublet(root->children[1], x, y, false);
    weight_t s0 = t0 == 0 ? 0 : aa_doublet(root->children[0], x, y);
    weight_t s1 = t1 == 0 ? 0 : aa_doublet(root->children[1], x, y);
    root->doublet[x][y] = s0 + s1 + root->children[0]->singlet[x] * t1 + root->children[1]->singlet[x] * t0;
    if (root->doublet[x][y] == 0) total_count[3] ++;
    return root->doublet[x][y];
}
*/

weight_t Tree::aa_doublet(Node *root, index_t x, index_t y) {
    if (root->children.size() == 0) return 0;
    weight_t c0 = root->children[0]->singlet[x], c1 = root->children[1]->singlet[x], s0 = 0, s1 = 0;
    weight_t t0 = get_doublet(root->children[0], x, y, false), t1 = get_doublet(root->children[1], x, y, false);
    if (c0 != 0 && t0 != 0) s0 = aa_doublet(root->children[0], x, y);
    if (c1 != 0 && t1 != 0) s1 = aa_doublet(root->children[1], x, y);
    weight_t c = s0 + s1 + c0 * t1 + c1 * t0;
    root->add_doublet(x, y, c);
    return c;
}

std::unordered_set<index_t> Tree::bad_edges(Node *root, Taxa &subset, weight_t ***graph) {
    if (root->children.size() == 0) {
        std::unordered_set<index_t> subtree;
        index_t index = subset.get_index(root->index);
        subtree.insert(index);
        return subtree;
    }
    else {
        std::unordered_set<index_t> l = bad_edges(root->children[0], subset, graph);
        std::unordered_set<index_t> r = bad_edges(root->children[1], subset, graph);
        std::unordered_set<index_t> subtree;
        for (auto i = l.begin(); i != l.end(); i ++) {
            if (subtree.find(*i) == subtree.end())
                subtree.insert(*i);
        }
        for (auto j = r.begin(); j != r.end(); j ++) {
            if (subtree.find(*j) == subtree.end()) 
                subtree.insert(*j);
        }
        for (auto i = l.begin(); i != l.end(); i ++) {
            for (auto j = r.begin(); j != r.end(); j ++) {
                if (*i == *j) continue;
                index_t x = subset.root_key(*i), y = subset.root_key(*j);
                index_t i_ = subset.root_index(*i), j_ = subset.root_index(*j);
                weight_t s = 0;
                if (x == 0) {
                    if (y == 0) {
                        s += index2node[*i]->get_doublet(0, 0) + index2node[*j]->get_doublet(0, 0);
                        s -= root->children[0]->get_doublet(0, 0) + root->children[1]->get_doublet(0, 0);
                        s += get_doublet(root, 0, 0, true);
                    }
                    else {
                        s += root->children[1]->get_doublet(y, 0);
                        weight_t t = root->children[1]->singlet[y];
                        s += t * (index2node[*i]->get_doublet(0, y) - root->children[0]->get_doublet(0, y));
                        s += t * get_doublet(root, y, 0, true);
                    }
                }
                else {
                    if (y == 0) {
                        s += root->children[0]->get_doublet(x, 0);
                        weight_t t = root->children[0]->singlet[x];
                        s += t * (index2node[*j]->get_doublet(0, x) - root->children[1]->get_doublet(0, x));
                        s += t * get_doublet(root, x, 0, true);
                    }
                    else {
                        weight_t t0 = root->children[0]->singlet[x], t1 = root->children[1]->singlet[y];
                        s += t1 * root->children[0]->get_doublet(x, y) + t0 * root->children[1]->get_doublet(y, x);
                        s += t1 * t0 * get_doublet(root, x, y, true);
                    }
                }
                graph[1][i_][j_] += s;
                graph[1][j_][i_] += s;
            }
        }
        return subtree;
    }
}

void Tree::good_edges(Node *root, Taxa &subset, weight_t ***graph) {
    index_t s = subset.singleton_taxa(), m = subset.artificial_taxa();
    std::unordered_map<index_t, index_t> tree_indices = get_indices();
    std::unordered_set<index_t> valid;
    for (auto elem : tree_indices) 
        valid.insert(subset.root_index(elem.first));
    if (subset.normalization() != '0') {
        weight_t *c = new weight_t[m + 1];
        for (index_t i = 0; i <= m; i ++) 
            c[i] = root->singlet[i];
        c[0] -= 2;
        weight_t sum = Node::get_doublet(c, root->s1, root->s2, 0, 0);
        delete [] c;
        for (index_t i = 0; i < subset.size(); i ++) {
            if (valid.find(i) == valid.end()) continue;
            for (index_t j = 0; j < subset.size(); j ++) {
                if (i == j || valid.find(j) == valid.end()) continue;
                graph[0][i][j] = sum - graph[1][i][j];
            }
        }
    }
    else {
        weight_t *c = new weight_t[m + 1];
        for (index_t i = 0; i <= m; i ++) 
            c[i] = root->singlet[i];
        c[0] -= 2;
        weight_t sum = Node::get_doublet(c, root->s1, root->s2, 0, 0);
        for (index_t i = 0; i < s; i ++) {
            if (valid.find(i) == valid.end()) continue;
            for (index_t j = 0; j < s; j ++) {
                if (i == j || valid.find(j) == valid.end()) continue;
                graph[0][i][j] = sum - graph[1][i][j];
            }
        }
        c[0] += 1;
        for (index_t i = 0; i < s; i ++) {
            if (valid.find(i) == valid.end()) continue;
            for (index_t j = 0; j < m; j ++) {
                if (valid.find(s + j) == valid.end()) continue;
                sum = Node::get_doublet(c, root->s1, root->s2, j + 1, 0) * root->singlet[j + 1];
                graph[0][s + j][i] = graph[0][i][s + j] = sum - graph[1][i][s + j];
            }
        }
        c[0] += 1;
        for (index_t i = 0; i < m; i ++) {
            if (valid.find(s + i) == valid.end()) continue;
            for (index_t j = 0; j < m; j ++) {
                if (i == j || valid.find(s + j) == valid.end()) continue;
                sum = Node::get_doublet(c, root->s1, root->s2, i + 1, j + 1);
                sum *= root->singlet[i + 1] * root->singlet[j + 1];
                graph[0][s + j][s + i] = graph[0][s + i][s + j] = sum - graph[1][s + i][s + j];
            }
        }
        delete [] c;
    }
}
