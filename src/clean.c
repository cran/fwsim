#include "clean.h"

static void clear_tree(struct kdtree *tree)
{
  struct kdnode* head;
  struct kdnode* current;
  struct kdnode* next;
  
  if (!tree || !tree->root) {
    return;
  }

  head = tree->root;
  
  if (head->right) {
    current = head->right;
      
    while (current != NULL) {
      next = current->right;
      free(current->pos);
      free(current);
      current = next;
      if (current == head) break;
    }
  }
  
	if (tree->rect) {
	  free(tree->rect->min);
    free(tree->rect->max);
	  free(tree->rect);
		tree->rect = 0;
	}
}

void free_tree(struct kdtree *tree)
{
	if (tree) {
		clear_tree(tree);
		free(tree);
	}
}

