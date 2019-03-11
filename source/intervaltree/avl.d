module intervaltree.avl;

/// https://rosettacode.org/wiki/AVL_tree#D

private struct IntervalTreeNode(IntervalType)
{
    IntervalType interval;

    // sort key
    alias key = interval.start;
    int balance;// balance factor
    int height; // aka size; # elements in subtree

    IntervalTreeNode *parent;
    IntervalTreeNode *left;
    IntervalTreeNode *right;

    invariant
    {
        // Ensure children are distinct
        if (this.left !is null && this.right !is null)
        {
            assert(this.left != this.right, "Left and righ child appear identical");
        }
    }
}

///
class IntervalAVLtree(IntervalType)
{
    alias Node = IntervalTreeNode!IntervalType;

    private Node* root;

    /**
    * Find a node in the tree
    *
    * @param x       node value to find (in)
    * @param cnt     number of nodes smaller than or equal to _x_; can be NULL (out)
    *
    * @return node equal to _x_ if present, or NULL if absent
    */
    Node* find(const(Node)* x, out uint cnt)
    {
        const Node* p = this.root;
        //uint cnt = 0;
        while (p !is null) {
            const int cmp = (x < p);
            if (cmp >= 0) cnt += (p.left ? p.left.height : 0) + 1;
            if (cmp < 0) p = p.left;
            else if (cmp > 0) p = p.right;
            else break;
        }
        //if (cnt_ !is null) *cnt_ = cnt;
        return p;
    }

    /// /* one rotation: (a,(b,c)q)p => ((a,b)p,c)q */ \
    /// /* dir=0 to left; dir=1 to right */
    private Node* rotate1(Node* p, const int dir)
    {
        const int opp = 1 - dir; /* opposite direction */
        Node *q = p.__head.p[opp];
        uint size_p = p.__head.size;
        p.__head.size -= q.__head.size - kavl_size_child(__head, q, dir);
        q.__head.size = size_p;
        p.__head.p[opp] = q.__head.p[dir];
        q.__head.p[dir] = p;
        return q;
    }

    ///
    pure nothrow @nogc @safe
    final bool insert(const Node key);

    ///
    pure nothrow @nogc @safe
    final bool remove(Node key);

    pure nothrow @nogc @safe
    private void rebalance(Node *n);


}


/* The MIT License
   Copyright (c) 2018 by Attractive Chaos <attractor@live.co.uk>
   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:
   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* An example:
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "kavl.h"
struct my_node {
  char key;
  KAVL_HEAD(struct my_node) head;
};
#define my_cmp(p, q) (((q)->key < (p)->key) - ((p)->key < (q)->key))
KAVL_INIT(my, struct my_node, head, my_cmp)
int main(void) {
  const char *str = "MNOLKQOPHIA"; // from wiki, except a duplicate
  struct my_node *root = 0;
  int i, l = strlen(str);
  for (i = 0; i < l; ++i) {        // insert in the input order
    struct my_node *q, *p = malloc(sizeof(*p));
    p->key = str[i];
    q = kavl_insert(my, &root, p, 0);
    if (p != q) free(p);           // if already present, free
  }
  kavl_itr_t(my) itr;
  kavl_itr_first(my, root, &itr);  // place at first
  do {                             // traverse
    const struct my_node *p = kavl_at(&itr);
    putchar(p->key);
    free((void*)p);                // free node
  } while (kavl_itr_next(my, &itr));
  putchar('\n');
  return 0;
}
*/

enum DIR : int
{
    LEFT = 0,
    RIGHT = 1
}

///
enum KAVL_MAX_DEPTH = 64;

///
pragma(inline, true)
@safe @nogc nothrow
auto kavl_size(T)(T* p) { return (p ? p.size : 0); }

///
pragma(inline, true)
@safe @nogc nothrow
auto kavl_size_child(T* q, DIR i) { return (q.p[i] ? q.p[i].size : 0); }

///
struct Node(T)
{
    T key;  /// sortable

    /// below prev in ->head
    Node*[2] p;
    char balance;   /// balance factor
    uint size;      /// #elements in subtree
}

/**
 * Find a node in the tree
 *
 * @param root    root of the tree
 * @param x       node value to find (in)
 * @param cnt     number of nodes smaller than or equal to _x_; can be NULL (out)
 *
 * @return node equal to _x_ if present, or NULL if absent
 */
Node *kavl_find(const(Node) *root, const(Node) *x, out uint cnt) {

    const(Node)* p = root;
    
    while (p !is null) {
        const int cmp = (x < p);
        if (cmp >= 0) cnt += kavl_size_child(p, DIR.LEFT) + 1; // left tree plus self

        if (cmp < 0) p = p.p[DIR.LEFT];         // descend leftward
        else if (cmp > 0) p = p.p[DIR.RIGHT];   // descend rightward
        else break;
    }

    return p;
}


/// /* one rotation: (a,(b,c)q)p => ((a,b)p,c)q */
pragma(inline, true)
@safe @nogc nothrow
private
Node *kavl_rotate1(Node *p, DIR dir) { /* dir=0 to left; dir=1 to right */
    const int opp = 1 - dir; /* opposite direction */
    Node *q = p.p[opp];
    const uint size_p = p.size;
    p.size -= q.size - kavl_size_child(q, dir);
    q.size = size_p;
    p.p[opp] = q.p[dir];
    q.p[dir] = p;
    return q;
}

/** two consecutive rotations: (a,((b,c)r,d)q)p => ((a,b)p,(c,d)q)r */
pragma(inline, true)
@safe @nogc nothrow
private
Node *kavl_rotate2(Node *p, DIR dir) {
    int b1;
    const int opp = 1 - dir;
    Node *q = p.p[opp], *r = q.p[dir];
    const uint size_x_dir = kavl_size_child(r, dir);
    r.size = p.size;
    p.size -= q.size - size_x_dir;
    q.size -= size_x_dir + 1;
    p.p[opp] = r.p[dir];
    r.p[dir] = p;
    q.p[dir] = r.p[opp];
    r.p[opp] = q;
    b1 = dir == 0 ? +1 : -1;
    if (r.balance == b1) q.balance = 0, p.balance = -b1;
    else if (r.balance == 0) q.balance = p.balance = 0;
    else q.balance = b1, p.balance = 0;
    r.balance = 0;
    return r;
}

/**
 * Insert a node to the tree
 *
 * @param suf     name suffix used in KAVL_INIT()
 * @param proot   pointer to the root of the tree (in/out: root may change)
 * @param x       node to insert (in)
 * @param cnt     number of nodes smaller than or equal to _x_; can be NULL (out)
 *
 * @return _x_ if not present in the tree, or the node equal to x.
 */
//#define kavl_insert(suf, proot, x, cnt) kavl_insert_##suf(proot, x, cnt)

//#define __KAVL_INSERT(suf, __scope, __type, __head, __cmp)
@safe @nogc nothrow
Node *kavl_insert(Node **root_, Node *x, out uint cnt)
{
    
    char[KAVL_MAX_DEPTH] stack;
    Node*[KAVL_MAX_DEPTH] path;

    Node* bp;
    Node* bq;
    Node* p;
    Node* q;
    Node* r = null; /* _r_ is potentially the new root */

    int i, which = 0, top, b1, path_len;

    bp = *root_, bq = 0;
    /* find the insertion location */
    for (p = bp, q = bq, top = path_len = 0; p; q = p, p = p.p[which]) {
        const int cmp = (x < p);
        if (cmp >= 0) cnt += kavl_size_child(p, 0) + 1;
        if (cmp == 0) {
            return p;
        }
        if (p.balance != 0)
            bq = q, bp = p, top = 0;
        stack[top++] = which = (cmp > 0);
        path[path_len++] = p;
    }

    x.balance = 0, x.size = 1, x.p[0] = x.p[1] = 0;
    if (q == 0) *root_ = x;
    else q.p[which] = x;
    if (bp == 0) return x;
    for (i = 0; i < path_len; ++i) ++path[i].size;
    for (p = bp, top = 0; p != x; p = p.p[stack[top]], ++top) /* update balance factors */
        if (stack[top] == 0) --p.balance;
        else ++p.balance;
    if (bp.balance > -2 && bp.balance < 2) return x; /* no re-balance needed */
    /* re-balance */
    which = (bp.balance < 0);
    b1 = which == 0 ? +1 : -1;
    q = bp.p[1 - which];
    if (q.balance == b1) {
        r = kavl_rotate1(bp, which);
        q.balance = bp.balance = 0;
    } else r = kavl_rotate2(bp, which);
    if (bq == 0) *root_ = r;
    else bq.p[bp != bq.p[0]] = r;
    return x;
}


/**
 * Delete a node from the tree
 *
 * @param suf     name suffix used in KAVL_INIT()
 * @param proot   pointer to the root of the tree (in/out: root may change)
 * @param x       node value to delete; if NULL, delete the first node (in)
 *
 * @return node removed from the tree if present, or NULL if absent
 */
/+
#define kavl_erase(suf, proot, x, cnt) kavl_erase_##suf(proot, x, cnt)
#define kavl_erase_first(suf, proot) kavl_erase_##suf(proot, 0, 0)

#define __KAVL_ERASE(suf, __scope, __type, __head, __cmp) \
	__scope __type *kavl_erase_##suf(__type **root_, const __type *x, unsigned *cnt_) { \
		__type *p, *path[KAVL_MAX_DEPTH], fake; \
		unsigned char dir[KAVL_MAX_DEPTH]; \
		int i, d = 0, cmp; \
		unsigned cnt = 0; \
		fake.__head.p[0] = *root_, fake.__head.p[1] = 0; \
		if (cnt_) *cnt_ = 0; \
		if (x) { \
			for (cmp = -1, p = &fake; cmp; cmp = __cmp(x, p)) { \
				int which = (cmp > 0); \
				if (cmp > 0) cnt += kavl_size_child(__head, p, 0) + 1; \
				dir[d] = which; \
				path[d++] = p; \
				p = p->__head.p[which]; \
				if (p == 0) { \
					if (cnt_) *cnt_ = 0; \
					return 0; \
				} \
			} \
			cnt += kavl_size_child(__head, p, 0) + 1; /* because p==x is not counted */ \
		} else { \
			for (p = &fake, cnt = 1; p; p = p->__head.p[0]) \
				dir[d] = 0, path[d++] = p; \
			p = path[--d]; \
		} \
		if (cnt_) *cnt_ = cnt; \
		for (i = 1; i < d; ++i) --path[i]->__head.size; \
		if (p->__head.p[1] == 0) { /* ((1,.)2,3)4 => (1,3)4; p=2 */ \
			path[d-1]->__head.p[dir[d-1]] = p->__head.p[0]; \
		} else { \
			__type *q = p->__head.p[1]; \
			if (q->__head.p[0] == 0) { /* ((1,2)3,4)5 => ((1)2,4)5; p=3 */ \
				q->__head.p[0] = p->__head.p[0]; \
				q->__head.balance = p->__head.balance; \
				path[d-1]->__head.p[dir[d-1]] = q; \
				path[d] = q, dir[d++] = 1; \
				q->__head.size = p->__head.size - 1; \
			} else { /* ((1,((.,2)3,4)5)6,7)8 => ((1,(2,4)5)3,7)8; p=6 */ \
				__type *r; \
				int e = d++; /* backup _d_ */\
				for (;;) { \
					dir[d] = 0; \
					path[d++] = q; \
					r = q->__head.p[0]; \
					if (r->__head.p[0] == 0) break; \
					q = r; \
				} \
				r->__head.p[0] = p->__head.p[0]; \
				q->__head.p[0] = r->__head.p[1]; \
				r->__head.p[1] = p->__head.p[1]; \
				r->__head.balance = p->__head.balance; \
				path[e-1]->__head.p[dir[e-1]] = r; \
				path[e] = r, dir[e] = 1; \
				for (i = e + 1; i < d; ++i) --path[i]->__head.size; \
				r->__head.size = p->__head.size - 1; \
			} \
		} \
		while (--d > 0) { \
			__type *q = path[d]; \
			int which, other, b1 = 1, b2 = 2; \
			which = dir[d], other = 1 - which; \
			if (which) b1 = -b1, b2 = -b2; \
			q->__head.balance += b1; \
			if (q->__head.balance == b1) break; \
			else if (q->__head.balance == b2) { \
				__type *r = q->__head.p[other]; \
				if (r->__head.balance == -b1) { \
					path[d-1]->__head.p[dir[d-1]] = kavl_rotate2_##suf(q, which); \
				} else { \
					path[d-1]->__head.p[dir[d-1]] = kavl_rotate1_##suf(q, which); \
					if (r->__head.balance == 0) { \
						r->__head.balance = -b1; \
						q->__head.balance = b1; \
						break; \
					} else r->__head.balance = q->__head.balance = 0; \
				} \
			} \
		} \
		*root_ = fake.__head.p[0]; \
		return p; \
	}
+/
/// free every node
void kavl_free(T)(auto __head, auto __root, auto __free) 
{
    do {
        T* _p;
        T* _q;
        for (_p = __root; _p; _p = _q) {
            if (_p.__head.p[0] == 0) {
                _q = _p.__head.p[1];
                __free(_p);
            } else {
                _q = _p.__head.p[0];
                _p.__head.p[0] = _q.__head.p[1];
                _q.__head.p[1] = _p;
            }
        }
    } while (0);
}
/+
#define __KAVL_ITR(suf, __scope, __type, __head, __cmp) \
	struct kavl_itr_##suf { \
		const __type *stack[KAVL_MAX_DEPTH], **top, *right; /* _right_ points to the right child of *top */ \
	}; \
	__scope void kavl_itr_first_##suf(const __type *root, struct kavl_itr_##suf *itr) { \
		const __type *p; \
		for (itr->top = itr->stack - 1, p = root; p; p = p->__head.p[0]) \
			*++itr->top = p; \
		itr->right = (*itr->top)->__head.p[1]; \
	} \
	__scope int kavl_itr_find_##suf(const __type *root, const __type *x, struct kavl_itr_##suf *itr) { \
		const __type *p = root; \
		itr->top = itr->stack - 1; \
		while (p != 0) { \
			int cmp; \
			cmp = __cmp(x, p); \
			if (cmp < 0) *++itr->top = p, p = p->__head.p[0]; \
			else if (cmp > 0) p = p->__head.p[1]; \
			else break; \
		} \
		if (p) { \
			*++itr->top = p; \
			itr->right = p->__head.p[1]; \
			return 1; \
		} else if (itr->top >= itr->stack) { \
			itr->right = (*itr->top)->__head.p[1]; \
			return 0; \
		} else return 0; \
	} \

/**
 * Place the iterator at the smallest object
 *
 * @param suf     name suffix used in KAVL_INIT()
 * @param root    root of the tree
 * @param itr     iterator
 */
#define kavl_itr_first(suf, root, itr) kavl_itr_first_##suf(root, itr)

/**
 * Place the iterator at the object equal to or greater than the query
 *
 * @param suf     name suffix used in KAVL_INIT()
 * @param root    root of the tree
 * @param x       query (in)
 * @param itr     iterator (out)
 *
 * @return 1 if find; 0 otherwise. kavl_at(itr) is NULL if and only if query is
 *         larger than all objects in the tree
 */
#define kavl_itr_find(suf, root, x, itr) kavl_itr_find_##suf(root, x, itr)
+/

/**
 * Move to the next object in order
 *
 * @param itr     iterator (modified)
 *
 * @return 1 if there is a next object; 0 otherwise
 */
int kavl_itr_next(T)(kavl_itr_t *itr) {
    for (;;) {
        const T *p;
        for (p = itr.right, --itr.top; p; p = p.__head.p[0])
            *++itr.top = p;
        if (itr.top < itr.stack) return 0;
        itr.right = (*itr.top).__head.p[1];
        return 1;
    }
}

/**
 * Return the pointer at the iterator
 *
 * @param itr     iterator
 *
 * @return pointer if present; NULL otherwise
 */
///#define kavl_at(itr) ((itr)->top < (itr)->stack? 0 : *(itr)->top)
pragma(inline, true)
auto kavl_at(kavl_itr_t *itr)
{
    return (itr.top < itr.stack) ? null : itr.top;
}