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
/* 
  Library modified heavily by Mikkel Meyer Andersen <mikl@math.aau.dk>
  for fwsim R package.
  A lot of clean up and changed applied. Also added
    haplotypes_equal
    kd_find_pos
    kd_insert_or_update_count
*/

#ifndef _KDTREE_H_
#define _KDTREE_H_

#ifdef __cplusplus
extern "C" {
#endif

struct kdhyperrect {
  int dim;
  int *min, *max;              /* minimum/maximum coords */
};

struct kdnode {
  int *pos;
  int dir;
  int count;

  struct kdnode *left, *right;  /* negative/positive side */
};

struct kdtree {
  int dim;
  struct kdnode *root;
  struct kdhyperrect *rect;
  void (*destr)(void*);
};

/* create a kd-tree for "k"-dimensional data */
struct kdtree *kd_create(int k);

/* insert a node, specifying its position and count */
int kd_insert(struct kdtree *tree, const int *pos, const int count);

/* insert a node or update existing */
int kd_insert_or_update_count(struct kdtree *tree, const int *pos, const int count);

#ifdef __cplusplus
}
#endif

#endif  /* _KDTREE_H_ */

