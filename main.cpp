#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <bits/stdc++.h> //for using log2 function
#include <map>
 

struct TreeNode;
struct EdgeNode;
 
typedef std::string tree_t;
 
struct EdgeNode{
    tree_t val;
    TreeNode* subtree;
    EdgeNode* next;
};
 
struct TreeNode{
    tree_t val;
    EdgeNode* subtree_l;
};

//it keeps only the rows specified by the string and column number parced into the function.
std::vector<std::vector<std::string>> keepRowsWithValue(const std::vector<std::vector<std::string>>& matrix, int columnIdx, const std::string& value) {
    std::vector<std::vector<std::string>> modifiedMatrix;
    // Keep the first row
    if (!matrix.empty()) {
        modifiedMatrix.push_back(matrix[0]);
    }
    for (size_t i = 1; i < matrix.size(); ++i) {
        const auto& row = matrix[i];
        if (row[columnIdx] == value) {
            modifiedMatrix.push_back(row);
        }
    }
    return modifiedMatrix;
}


//swaps the rows of a 2D matrix 
void swap_rows(std::vector<std::vector<std::string>> &matrix, int row1, int row2){
    std::vector<std::string> temp = matrix[row1];
    matrix[row1] = matrix[row2];
    matrix[row2] = temp;
}


//bubble sort of rows in matrix depending on column chosen. Without break statement which stops sorting once no changes are detected.
void bubble_sort(std::vector<std::vector<std::string>> &matrix, int column){
    int msize = matrix.size();
    for(int i = 0; i <msize-1; i++){ //total number of iterations for bubble sort
        for(int j = 1; j < msize-1-i; j++){ //columns, sorts for 1 less element each loop as the last element after each loop is always largest
            if(matrix[j][column] > matrix[j+1][column]){
                swap_rows(matrix, j, j+1);
            }
        }
    }
}


//counts the total number of different options available as an awser for each attirbute, and then lists how many time that option has been represented
//stores values in 2D array
//ie: tempreture: 4 choices are low, 6 choices are high. value stored as [4,6] in terms of alphabetical order in 2D array.
std::vector<std::vector<int>> count_inputs(std::vector<std::vector<std::string>> matrix){

    // Vector to store the counts of different input options for each attribute
    std::vector<std::vector<int>> v_attributes; 
    // Iterate over each row in the matrix
    for (int j = 0; j < matrix[0].size(); j++){
        int count_weight_of_anwsers = 1; //starting with 1 as first value. First row is title names
        bubble_sort(matrix, j); 

        //total number of times each different anwsers pops up for each attribute
        std::vector<int> v_amount; 

        //itterate over each row of matrix except 1.
        for (int i = 1; i < matrix.size(); i++){

            //if adjacent elements the same then increment count of that number
            if (i+1 < matrix.size() && matrix[i][j] == matrix[i+1][j]){
                count_weight_of_anwsers++;
            }

            // Push the count of the previous set of answers and reset the count if adjacent elements arent identical
            if (i+1 < matrix.size() && matrix[i][j] != matrix[i+1][j]){
                v_amount.push_back(count_weight_of_anwsers);
                count_weight_of_anwsers = 1;
            }
        }
        v_amount.push_back(count_weight_of_anwsers);
        v_attributes.push_back(v_amount);
    }
    return(v_attributes); 
}

//returns the possible choices of a column in a vector
//example: first column tempreture has onle two choices ={high, low}
std::vector<std::string> choices_of_column(std::vector<std::vector<std::string>> matrix, int column_n){
    std::vector<std::string> choices_vector;
    bubble_sort(matrix, column_n); //sort and add a new element every time adjacent elements arent the same
    for (int i = 0; i < matrix.size(); i++) {
        if (i+1 < matrix.size() && matrix[i][column_n] != matrix[i+1][column_n]){
            choices_vector.push_back(matrix[i+1][column_n]);
        }
    }
    return(choices_vector);
}

//counts the number of times each output occurs with each split
//saved in a 2D vector, each element representing a single outcome of the split 
//example: tempreture will have 2 maps representing it. 1 map for high and 1 map for low
//so temp vector will look like: tempreture = {{acceptable: 2, good: 2, poor: 2}, {acceptable: 1, good: 1, poor: 2}} 
//the first map is for high and second is for low
std::vector<std::map<std::string, int>> count_outcomes_after_split(std::vector<std::vector<std::string>> matrix, int column_n){
    std::vector<std::map<std::string, int>> attribute_vector; //stored outcome for each split
    std::map<std::string, int> attribute_map;//map storing specific choice: ie: high

    bubble_sort(matrix, column_n); 

    //itterate through the whole matrix to count the outcomes after split (itterate through every row of a chosen column)
    for (int i = 1; i < matrix.size(); i++) {

        // If the current value is the same as the next value in the split column
        if (i+1 < matrix.size() && matrix[i][column_n] == matrix[i+1][column_n]){
            auto it = attribute_map.find(matrix[i][matrix[0].size()-1]) ;
            if (it != attribute_map.end()) {
                // If the outcome already exists in the attribute map, increment its count
                attribute_map[matrix[i][matrix[0].size()-1]]++;
            }
            else {
                // If the outcome doesn't exist in the attribute map, add it with a count of 1
                attribute_map[matrix[i][matrix[0].size()-1]] = 1;
            }
        }

        // If the current value is different from the next value in the split column
        if (i+1 < matrix.size() && matrix[i][column_n] != matrix[i+1][column_n]){
            auto it = attribute_map.find(matrix[i][matrix[0].size()-1]) ;
            if (it != attribute_map.end()) {
                // If the outcome already exists in the attribute map, increment its count
                attribute_map[matrix[i][matrix[0].size()-1]]++;
            }
            else {
                // If the outcome doesn't exist in the attribute map, add it with a count of 1
                attribute_map[matrix[i][matrix[0].size()-1]] = 1;
            }

            // Add the attribute map to the attribute vector and clear the map for the next split
            attribute_vector.push_back(attribute_map);
            attribute_map.clear();
        }

        //if you have reached the last element in the matrix.
        if (i+1 == matrix.size()){
            auto it = attribute_map.find(matrix[i][matrix[0].size()-1]) ;
            if (it != attribute_map.end()) {
                // If the outcome already exists in the attribute map, increment its count
                attribute_map[matrix[i][matrix[0].size()-1]]++;
            }
            else {
                // If the outcome doesn't exist in the attribute map, add it with a count of 1
                attribute_map[matrix[i][matrix[0].size()-1]] = 1;
            }
            // Add the attribute map to the attribute vector and clear the map for the next split
            attribute_vector.push_back(attribute_map);
            attribute_map.clear();
        }

    }
    
    
    return (attribute_vector);
}

//finds the entropy of every single element in a column and returns in a vector of doubles.
//used to find shortcutes. ie: if any element in this vector is 0, we will know that it is completely deterministic.
//therefore instead of splitting fuether, we go straight to the output
//example: wind --(strong)--> poor, because entropy of poor is 0.
//not taking weight into account as it is not used for the pupose of finding next split
std::vector<double> entropy_of_column(std::vector<std::vector<std::string>> matrix, int column_n){

    //get all different choices 
    std::vector<std::vector<int>> c_inputs = count_inputs(matrix);
    std::vector<double> entropy_final; //information gain vector sotring all possible information gains after each split.
    double e_inputs = 0.0;

    // Vector of maps representing outcomes after each split
    std::vector<std::map<std::string, int>> c_outcomes = count_outcomes_after_split(matrix, column_n); 
    std::vector<int> input_occurances = c_inputs[column_n]; //map of how many occurances of inpts there are for split

    // Calculate entropy for each input element in the column
    for (int j = 0; j < input_occurances.size(); j++){
        e_inputs = 0.0;
        std::map<std::string, int> map = c_outcomes[j];

        //itterating over each element in the map
        for (const auto& pair : map) {
            int value = pair.second;

            // Calculate the entropy of the outcome occurrence
            double x = static_cast<double>(value) / input_occurances[j] ;
            // Calculate the entropy contribution and subtract it from the total entropy
            e_inputs -= x*log2(x);        
        }
        entropy_final.push_back(e_inputs);
    }
    return(entropy_final);

}

//returns in a vector the entropies of each column except the outcome column
//example: {0.654, 1.534, 0.876}
//used to pick which column to split by next
std::vector<double> entropy_of_matrix(std::vector<std::vector<std::string>> matrix){

    std::vector<std::vector<int>> c_inputs = count_inputs(matrix);
    std::vector<double> entropy_final; //information gain vector sotring all possible information gains after each split.
    int tot = matrix.size() - 1; //used for finding entropy 

    //calculate the entropy after each split, for each column
    for (int i = 0; i < matrix[0].size()-1; i++){

        std::vector<std::map<std::string, int>> c_outcomes = count_outcomes_after_split(matrix, i); //vector of maps for each column
        std::vector<int> input_occurances = c_inputs[i]; //map of how many occurances of inpts there are for split
        double e_inputs_tot = 0.0;

        for (int j = 0; j < input_occurances.size(); j++){

            double e_inputs = 0.0;
            std::map<std::string, int> map = c_outcomes[j];

            //itterating over each element in the map
            for (const auto& pair : map) {
                int value = pair.second;

                //calculate entropy of each outcome
                double p = static_cast<double>(value) / input_occurances[j] ;
                e_inputs -= p*log2(p);

            }

            //calculate the weight of each split and add to total entropy value for the split of that input
            double weight = static_cast<double>(input_occurances[j]) / (tot);

            e_inputs_tot +=  weight * e_inputs;

        }
        entropy_final.push_back(e_inputs_tot);
    }
    return(entropy_final);
}




class A3Tree{

public:
    TreeNode* return_root(){
        TreeNode* root = t;
        return (root);
    }

    std::vector<int> return_order(){
        std::vector<int> order1 = order;
        return (order1);
    }

    TreeNode* allocate_tree_node(tree_t e){
            TreeNode* tmp = new TreeNode;
            tmp->val = e;
            tmp->subtree_l = NULL;
            return tmp; 
        }

    EdgeNode* cons_edge_node(TreeNode* t, EdgeNode* subtree_l, tree_t v){
        EdgeNode* tmp = new EdgeNode;
        tmp->subtree = t;
        tmp->next = subtree_l;
        tmp->val = v;
        return tmp;
    }

    void printTree(TreeNode* node, int indent = 0) {
        if (node == NULL) {
            return;
        }

        // Print the current node with proper indentation
        for (int i = 0; i < indent; i++) {
            std::cout << "  ";
        }
        std::cout << node->val << std::endl;

        // Recursively print the subtree edges
        EdgeNode* edge = node->subtree_l;
        while (edge != NULL) {
            printTree(edge->subtree, indent + 1);
            edge = edge->next;
        } 
    }

    //counts number of leaf nodes
    int CountLeafNodes(TreeNode* node) {
        if (node == NULL) {
            return 0; //not leaf
        }
        if (node->subtree_l == NULL) {
            return 1; //leaf
        }

        //itterates over each branch
        int leafCount = 0;
        EdgeNode* edge = node->subtree_l;
        while (edge != NULL) {
            leafCount += CountLeafNodes(edge->subtree);
            edge = edge->next;
        }

        return leafCount;
    }

    //leaf node assistant to avoid data hiding while keeping private data unchanged
    int leaf_node_count() {
        return CountLeafNodes(return_root());
    }

    int CountNodes(TreeNode* node) {
        if (node == nullptr) {
            return 0;
        }
        int totalCount = 1;  // Count the current node
        EdgeNode* edge = node->subtree_l;
        while (edge != nullptr) {
            totalCount += CountNodes(edge->subtree);
            edge = edge->next;
        }
        return totalCount;
    }

    //node node assistant to avoid data hiding while keeping private data unchanged
    int node_count() {
        return CountNodes(return_root());
    }


    //function used to keep track of order of elements, later used in query
    void swapElements(std::vector<int>& vec, int index1, int index2) {
        // Check if the indices are within the valid range of the vector
        if (index1 >= 0 && index1 < vec.size() && index2 >= 0 && index2 < vec.size()) {
            //swap elements
            std::swap(vec[index1], vec[index2]);
        } else {
            std::cout << "Please input exisiting elements for swap." << std::endl;
        }
        index1++; //each swap the next element will be swapped around (first swap 3rd column was swapped with first, second 4th column swapped with second and so on.)
    }


    //used in query, swaps the order of the strings provided in accordance to the order of parent ndos in tree
    //example: v = {"1", "2", "3"} will become {"1", "3", "2"} if order = {1, 3, 2}
    std::vector<std::string> swap_string_order(std::vector<std::string> input, const std::vector<int>& order){

        //takes in as const as it is private member data so swappedString defined .
        std::vector<std::string> swappedString = input;

        //itterates exactly rounded up(order.size/2). Because if it itterated over the whole of order.size(), we would scramble then 
        //unscramble the matrix. to avoid we only do half of the matrix. at odd numbers we round up.
        for (int i = 0; i < std::ceil(static_cast<double>(order.size()) / 2); i++) {
            std::string temp = swappedString[i];
            swappedString[i] = swappedString[order[i]];
            swappedString[order[i]] = temp;
            
        }
        return swappedString;


    }

    //used inside query to query
    tree_t query_super(TreeNode* root, std::vector<std::string> input) {

        tree_t found; //value found
        EdgeNode* it = root -> subtree_l; //itterated over child node of current node
        TreeNode* newParent = root; //tracks new parent node

        //check if at an edge (output) value, if no, continue the recusrsion
        if (newParent ->subtree_l != NULL){
            for (int i = 0; i < input.size(); i++){
                while (it != NULL){ //while not reached end of edgenodes
                    if (it -> val == input[i]){
                        input.erase(input.begin()); //remove first elements from input vector to itterate over next vector.
                        newParent = it ->subtree; 
                        found = query_super(newParent, input);                        
                    }
                    it = it->next;
                }
                return found;
            }
        }

        //if current parent node doesnt have child node return its value.
        if (newParent->subtree_l == NULL){
            return newParent->val;
        }
        return newParent->val;
    } 

    tree_t query(std::vector<std::string> input) {
        std::vector<std::string> swappedString = swap_string_order(input, return_order());
        return query_super(return_root(), swappedString);
    }

    TreeNode* add_tree_node(TreeNode* parent_node, tree_t child_node, tree_t edge_node) {
        auto parent = parent_node;

        if (parent == nullptr) {
            return NULL;
        }
        auto child = allocate_tree_node(child_node);
        auto new_edge = cons_edge_node(child, nullptr, edge_node); // Initialize next pointer as nullptr

        // If parent node doesn't have any edge nodes yet
        if (parent->subtree_l == nullptr) {
            parent->subtree_l = new_edge;
        } else {
            // Find the last edge node in the list and append the new edge node
            EdgeNode* last_edge = parent->subtree_l;
            while (last_edge->next != nullptr) {
                last_edge = last_edge->next;
            }
            last_edge->next = new_edge;
        }
        return(child);
    }



    TreeNode* build_tree_root(tree_t e){
        return allocate_tree_node(e);
    }


    TreeNode* build_tree(std::vector<std::vector<std::string>> matrix, TreeNode* parent_node, tree_t edge_val){

        //finding which column and rows to remove for new matrix 
        int index;
        for (int i = 0; i < matrix[0].size(); i++){
            if (parent_node->val == matrix[0][i]){
                index = i;
            }
        }
        //newMatrix has only needed rows and columns for next level of calculation: ie: once wind has been placed to the top we run the same entropy calcs with a smaller matrix
        std::vector<std::vector<std::string>> newMatrix = keepRowsWithValue(matrix, index, edge_val); //removing columns and rows for new matrix

        //entropy of new matrix columns
        std::vector<double> entropy_matrix = entropy_of_matrix(newMatrix);
        
        int minIndex = 0;
        double minValue = entropy_matrix[0];
        //finding which column has lowest entropy out of all.
        for (int i = 1; i < entropy_matrix.size(); i++) {
            if (entropy_matrix[i] < minValue) {
                minValue = entropy_matrix[i];
                minIndex = i;
            }
        }

        swapElements(order, swap_number, minIndex); //records order

        TreeNode* t_parent;
        t_parent = add_tree_node(parent_node, matrix[0][minIndex], edge_val);
        std::vector<double> entropy_column_chosen = entropy_of_column(newMatrix, minIndex);

        for (int i = 0; i < entropy_column_chosen.size(); i++){

            std::vector<std::string> choices = choices_of_column(newMatrix, minIndex);
            std::vector<std::map<std::string, int>> map_splits = count_outcomes_after_split(newMatrix, minIndex);

            if (entropy_column_chosen[i] != 0){
                TreeNode* t_new = build_tree(newMatrix, t_parent, choices[i]);
            }

            if (entropy_column_chosen[i] == 0){
                std::map<std::string, int> map = map_splits[i];
                std::string key = map.begin()->first;
                TreeNode* a = add_tree_node(t_parent, key, choices[i]);
            }
        }

    }



    //consutrcutor
    A3Tree (std::vector<std::vector<std::string>> matrix){
        std::vector<double> entropy_matrix = entropy_of_matrix(matrix);

        //printing total number of elements in matrix[0] and storing as private member data. used for query later. 
        for (int i = 0; i < matrix[0].size()-1; i++ ){ 
            order.push_back(i);
        }

        int minIndex = 0;
        double minValue = entropy_matrix[0];
        //finding which column has lowest entropy out of all.
        for (int i = 1; i < matrix[0].size()-1; i++) {
            if (entropy_matrix[i] < minValue) {
                minValue = entropy_matrix[i];
                minIndex = i;
            }
        }

        swapElements(order, swap_number, minIndex); //records order

        //finding individual entropies of choices in column to see if any have entropy 0.
        std::vector<double> entropy_column_chosen = entropy_of_column(matrix, minIndex);

        //create parent node: on the first itteration of code this is the root node.
        TreeNode* t1;
        t1 = build_tree_root(matrix[0][minIndex]);
        t = t1;
        for (int i = 0; i < entropy_column_chosen.size(); i++){ //for number of different choices possible in a column

            std::vector<std::string> choices = choices_of_column(matrix, minIndex);
            std::vector<std::map<std::string, int>> map_splits = count_outcomes_after_split(matrix, minIndex);

            //if entropy of choice is != 0 (therefore we need to build another parent node)
            if (entropy_column_chosen[i] != 0.0){
               TreeNode* t_new = build_tree(matrix, t1, choices[i]);

            }

            //check if any choices have an entropy of 0. if so then just construct the outcome as a parent node.
            if (entropy_column_chosen[i] == 0){
                std::map<std::string, int> map = map_splits[i];
                std::string key = map.begin()->first;
                TreeNode* a = add_tree_node(t1, key, choices[i]);
            }
        }
    }

    //destructor
    void deallocate_tree(TreeNode* t){
        if (t == nullptr)
            return;

        auto it = t->subtree_l;
        while (it != nullptr)
        {
            deallocate_tree(it->subtree);
            auto prev = it;
            it = it->next;
            delete prev;
        }
        
        delete t;
    }
    
private:
 
    TreeNode* t;
    std::vector<int> order;
    int swap_number = 0; 
};


int main(){
 
    // direct initialisation of a vector
    // in this case it's a vector containing vectors
    // each of which contains words (std::string)
    std::vector<std::vector<std::string>> matrix
    {
        {"temperature", "rain", "wind", "quality"},
        {"high", "yes", "light", "acceptable"},
        {"low", "yes", "light", "acceptable"},
        {"low", "no", "moderate", "good"},
        {"high", "yes", "strong", "poor"},
        {"high", "yes", "moderate", "acceptable"},
        {"high", "no", "moderate", "good"},
        {"low", "yes", "strong", "poor"},
        {"high", "no", "light", "good"},
        {"low", "yes", "moderate", "poor"},
        {"high", "no", "strong", "poor"},

    };
 
    A3Tree t(matrix);
    t.printTree(t.return_root());
    
    std::vector<std::string> q =  {"high", "yes", "light"};
    std::cout << t.query(q) << std::endl;
     
    std::cout << t.node_count() << std::endl;

    std::cout << t.leaf_node_count() << std::endl;

    t.deallocate_tree(t.return_root());

}