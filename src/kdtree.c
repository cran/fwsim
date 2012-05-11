/*
This file is part of ``kdtree'', a library for working with kd-trees.
Copyright (C) 2007-2011 John Tsiombikas <nuclear@member.fsf.org>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/
/* single nearest neighbor search written by Tamas Nepusz <tamas@cs.rhul.ac.uk> */
/* 
  Library modified heavily by Mikkel Meyer Andersen <mikl@math.aau.dk>
  for fwsim R package.
  A lot of clean up and changed applied. Also added
    haplotypes_equal
    kd_find_pos
    kd_insert_or_update_count
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "common.h"
#include "kdtree.h"
#include "print.h"

#if defined(WIN32) || defined(__WIN32__)
#include <malloc.h>
#endif

#define ABS(x)           (((x) < 0) ? -(x) : (x))

static int insert_rec(struct kdnode **node, const int *pos, const int count, int dir, int dim);

static struct kdhyperrect* hyperrect_create(int dim, const int *min, const int *max);
static void hyperrect_free(struct kdhyperrect *rect);
static struct kdhyperrect* hyperrect_duplicate(const struct kdhyperrect *rect);
static void hyperrect_extend(struct kdhyperrect *rect, const int *pos);
static int hyperrect_dist(struct kdhyperrect *rect, const int *pos);

#define alloc_resnode()    malloc(sizeof(struct res_node))
#define free_resnode(n)    free(n)

struct kdtree *kd_create(int k)
{
  struct kdtree *tree;

  if(!(tree = malloc(sizeof *tree))) {
    return 0;
  }

  tree->dim = k;
  tree->root = 0;
  tree->rect = 0;
  tree->size = 0;

  return tree;
}

static int insert_rec(struct kdnode **nptr, const int *pos, const int count, int dir, int dim)
{
  int new_dir;
  struct kdnode *node;

  if(!*nptr) {
    
    if(!(node = malloc(sizeof *node))) {
      return -1;
    }
    if(!(node->pos = malloc(dim * sizeof *node->pos))) {
      free(node);
      return -1;
    }
    memcpy(node->pos, pos, dim * sizeof *node->pos);

    node->count = count;
    node->dir = dir;
    node->left = node->right = 0;
    *nptr = node;
    return 0;   
  }

  node = *nptr;
  new_dir = (node->dir + 1) % dim;
  if(pos[node->dir] < node->pos[node->dir]) {
    return insert_rec(&(*nptr)->left, pos, count, new_dir, dim);
  }
  return insert_rec(&(*nptr)->right, pos, count, new_dir, dim);
}

int kd_insert(struct kdtree *tree, const int *pos, const int count)
{ 
  if (insert_rec(&tree->root, pos, count, 0, tree->dim)) {
    return -1;
  }

  tree->size += 1;
      
  if (tree->rect == 0) {
    tree->rect = hyperrect_create(tree->dim, pos, pos);
  } else {
    hyperrect_extend(tree->rect, pos);
  }

  return 0;
}

static int haplotypes_equal(const int dim, const int* h1, const int* h2) {
  int i;
  
  for (i = 0; i < dim; i++) {
    if (h1[i] != h2[i]) {
      return 0;
    }
  }
  
  return 1;
}

static struct kdnode* kd_find_pos(struct kdnode *node, const int *pos, int dir, int dim)
{ 
  int new_dir;
  
  if (!node) {
    return NULL;
  }
  
  /*
  printf("Comparing ");
  print_h(pos, dim);
  printf(" with ");
  print_h(node->pos, dim);
  printf("\n");
  */
  
  if (haplotypes_equal(dim, pos, node->pos)) {
    return node;
  }

  new_dir = (node->dir + 1) % dim;
  
  if(pos[node->dir] < node->pos[node->dir]) {
    return kd_find_pos(node->left, pos, new_dir, dim);
  }
  
  return kd_find_pos(node->right, pos, new_dir, dim);
}

int kd_insert_or_update_count(struct kdtree *tree, const int *pos, int count)
{
  struct kdres *res;
  int* e_haplotype;
  int e_is_equal;
  int e_dist;
 
  struct kdnode* node = tree->root;
  
  if (count <= 0) {
    error("Got count = 0, wasn't expected to reach such point...\n");
  }
  
  node = kd_find_pos(node, pos, 0, tree->dim);

  if (node) {
    node->count += count;
    return 0;
  }
  
  return kd_insert(tree, pos, count);
}

/* ---- hyperrectangle helpers ---- */
static struct kdhyperrect* hyperrect_create(int dim, const int *min, const int *max)
{
  size_t size = dim * sizeof(int);
  struct kdhyperrect* rect = 0;

  if (!(rect = malloc(sizeof(struct kdhyperrect)))) {
    return 0;
  }

  rect->dim = dim;
  if (!(rect->min = malloc(size))) {
    free(rect);
    return 0;
  }
  if (!(rect->max = malloc(size))) {
    free(rect->min);
    free(rect);
    return 0;
  }
  memcpy(rect->min, min, size);
  memcpy(rect->max, max, size);

  return rect;
}

static void hyperrect_free(struct kdhyperrect *rect)
{
  free(rect->min);
  free(rect->max);
  free(rect);
}

static struct kdhyperrect* hyperrect_duplicate(const struct kdhyperrect *rect)
{
  return hyperrect_create(rect->dim, rect->min, rect->max);
}

static void hyperrect_extend(struct kdhyperrect *rect, const int *pos)
{
  int i;

  for (i=0; i < rect->dim; i++) {
    if (pos[i] < rect->min[i]) {
      rect->min[i] = pos[i];
    }
    if (pos[i] > rect->max[i]) {
      rect->max[i] = pos[i];
    }
  }
}

static int hyperrect_dist(struct kdhyperrect *rect, const int *pos)
{
  int i;
  int result = 0;

  for (i=0; i < rect->dim; i++) {
    if (pos[i] < rect->min[i]) {
      result += ABS(rect->min[i] - pos[i]);
    } else if (pos[i] > rect->max[i]) {
      result += ABS(rect->max[i] - pos[i]);
    }
  }

  return result;
}

