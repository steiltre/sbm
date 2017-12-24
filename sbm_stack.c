
#include <stdlib.h>
#include "sbm_stack.h"

sbm_ll_node * sbm_ll_node_alloc()
{
  sbm_ll_node * new_node = malloc(sizeof*(new_node));

  new_node->val = -1;
  new_node->next = NULL;

  return new_node;
}

void sbm_ll_node_push(
    int32_t val,
    sbm_ll_node * root)
{
  sbm_ll_node * new_node = sbm_ll_node_alloc();

  new_node->val = val;
  new_node->next = root->next;
  root->next = new_node;
}

int32_t sbm_ll_node_pop(
    sbm_ll_node * root)
{
  sbm_ll_node * top = root->next;
  root->next = top->next;
  int32_t ret = top->val;
  free(top);
  return ret;
}
