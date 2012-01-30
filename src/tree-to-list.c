#include "tree-to-list.h"

/*
http://cslibrary.stanford.edu/109/TreeListRecursion.html
*/
static void join(struct kdnode *a, struct kdnode *b) {
    a->right = b;
    b->left = a;
}

static struct kdnode* append(struct kdnode *a, struct kdnode *b) {
    struct kdnode *aLast;
    struct kdnode *bLast;
    
    if (a == NULL) return(b);
    if (b == NULL) return(a);
    
    aLast = a->left;
    bLast = b->left;
    
    join(aLast, b);
    join(bLast, a);
    
    return(a);
}

/*
 --Recursion--
 Given an ordered binary tree, recursively change it into
 a circular doubly linked list which is returned.
*/
static struct kdnode* list_from_root_node(struct kdnode *root) {
    struct kdnode *aList;
    struct kdnode *bList;
    
    if (root == NULL) return(NULL);

    /* recursively solve subtrees -- leap of faith! */
    aList = list_from_root_node(root->left);
    bList = list_from_root_node(root->right);
    
    /* Make a length-1 list ouf of the root */
    root->left = root;
    root->right = root;

    /* Append everything together in sorted order */
    aList = append(aList, root);
    aList = append(aList, bList);
    
    return(aList);
}

struct kdnode* tree_to_list(const struct kdtree *tree, int r) {
  struct kdnode *head;  
  head = list_from_root_node(tree->root);
  return head;
}

