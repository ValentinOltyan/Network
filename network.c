#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <float.h>

const int MAXI = 11;
const int MAXD = 310;
struct pair{
  int in, out;
};
typedef struct pair pair;

struct node {
  int value;
  double w;
  struct node *next;
};
typedef struct node node;

struct graph {
  int size, r_size, nsize, nr_size;
  int (*nodes)[];
  node (*gr)[];
};
typedef struct graph graph;

const node empty = (node){0,0,NULL};
//create a graph
graph *newGraph(){
  graph *g = malloc(sizeof(graph));
  node (*n)[] = malloc(24 * sizeof(node));
  int (*a)[] = malloc(24 * sizeof(int));
  for(int i = 0; i < 24; i++){
    (*n)[i] = empty;
    (*a)[i] = -1;
  }
  g->nodes = a;
  g->gr = n;
  g->size = 24;
  g->r_size = 0;
  g->nsize = 24;
  g->nr_size = 0;
  return g;
}
//create a node
node *newNode(int val, double w){
  node *n = malloc(sizeof(node));
  n->value = val;
  n->w = w;
  n->next = NULL;
  return n;
}

bool isEmpty(node a){
  if(a.value == 0 && a.w == 0 && a.next == NULL)
    return true;
  return false;
}

bool isNode(graph *g, int x){
  for(int i = 0; i < g->nr_size; i++){
    if(x == (*g->nodes)[i])
      return true;
  }
  return false;
}

//pushes an edge in the graph
void push(graph* g, int poz, int val, double w){
  while(g->size <= g->r_size){
    g->gr = realloc(g->gr, (g->size * 3 / 2) * sizeof(node));
    g->size = g->size * 3 / 2;
  }
  if(!isNode(g, poz)){
    while(g->nsize <= g->nr_size){
      g->nodes = realloc(g->nodes, g->nsize * 3 / 2 * sizeof(int));
      g->nsize = g->nsize * 3 / 2;
    }
    (*g->nodes)[g->nr_size] = poz;
    g->nr_size ++;
    (*g->gr)[poz] = empty;
  }
  if(!isNode(g, val)){
    while(g->nsize <= g->nr_size){
      g->nodes = realloc(g->nodes, g->nsize * 3 / 2 * sizeof(int));
      g->nsize = g->nsize * 3 / 2;
    }
    (*g->nodes)[g->nr_size] = val;
    g->nr_size ++;
    (*g->gr)[val] = empty;
  }

  if(isEmpty((*g->gr)[poz])){
    (*g->gr)[poz] = (node){val,w,NULL};
    return;
  }
  node *a;
  a = &(*g->gr)[poz];
  while(a->next != NULL){
    a = a->next;
  }
  a->next = newNode(val, w);
}

// compares two strings that represent integers
int compare(const char a[], const char b[]){
  long n = strlen(a);
  long m = strlen(b);
  if(n > m)
    return 1;
  else
    if(n < m)
      return -1;
    for(long i = 0; i < n; i++){
        if(a[i] < b[i])
          return -1;
          else
          if(a[i] > b[i])
            return 1;
    }
  return 0;
}

//verifies if the given string is less than 50000000
bool verifStr(char s[]){
  if(compare(s,"50000000") > 0)
    return false;
  return true;
}
//verifies if the given string representing a double is too big
bool verifStrW(char s[]){
  if(strlen(s) > 10)
    return false;
  return true;
}
//handles the arguments given with the program call
void argDec(graph *g, char s[], int *poz, int *val, double *w){
  char a1[MAXI], a2[MAXI], a3[MAXD];
  if(!strchr(s,'-')){
    int a = atoi(s);
    if(!isNode(g, a)){
      if(g->r_size <= a)
        g->r_size = a;
      while(g->nsize <= a){
        g->nodes = realloc(g->nodes, (g->nsize * 3 / 2) * sizeof(int));
        g->nsize = g->nsize * 3 / 2;
      }
      while(g->size <= a){
        g->gr = realloc((*g->gr), (g->size * 3 / 2) * sizeof(node));
        g->size = g->size * 3 / 2;
      }
      (*g->nodes)[g->nr_size] = a;
      g->nr_size ++;
      (*g->gr)[a] = empty;
    }
  }
  else{
    sscanf(s,"%[^-]-%[^/]/%s",a1,a2,a3);
    if(verifStr(a1) != true){
      printf("Invalid node value \n");
      exit(1);
    }
    if(verifStr(a2) != true){
      printf("Invalid node value \n");
      exit(1);
    }
    if(strlen(a1) + strlen(a2) + 1 == strlen(s))
      *w = 1;
    else{
      if(verifStrW(a3) != true){
        printf("Invalid weight value \n");
        exit(1);
      }
      int k = 0;
      for(int i = 0 ; i< strlen(a3); i++)
      {
        if(!((a3[i] >= '0' && a3[i] <= '9') || (a3[i] == '.' && k == 0))){
          printf("Invalid input\n");
          exit(1);
        }
        if(a3[i] == '.')
          k = 1;
      }
      *w = atof(a3);
    }
    for(int i = 0 ; i< strlen(a1); i++)
    {
      if(!(a1[i] >= '0' && a1[i] <= '9')){
        printf("Invalid node value\n");
        exit(1);
      }
    }
    for(int i = 0 ; i< strlen(a2); i++)
    {
      if(!(a2[i] >= '0' && a2[i] <= '9')){
        printf("Invalid node value\n");
        exit(1);
      }
    }
    *poz = atoi(a1);
    *val = atoi(a2);
    if(*poz == 0 || *val == 0){
      printf("Node 0 is invalid\n");
      exit(1);
    }
  }
}
// comp function for qsort
int comp(const void *p, const void *f){
  const int *pi = p, *fi = f;
  int a = *pi, b = *fi;
  if(a > b)
    return 1;
  else
    if(a == b)
      return 0;
  return -1;
}
// create a network
graph *createNetwork(char s[]){
  char *p;
  graph *g = newGraph();
  p = strtok(s,",");
  do{
    int poz = 0, val = 0;
    double w = 0;
    argDec(g, p, &poz, &val, &w);
    if(g->r_size < poz)
      g->r_size = poz;
    if(g->r_size < val)
      g->r_size = val;
    if(poz != 0 && val != 0)
      push(g, poz, val, w);
    p = strtok(NULL,",");
  }while(p);
  qsort(*g->nodes, g->nr_size, sizeof(int), comp);
  return g;
}

//free a node and all its neighbours
void freeNode(node *a){
    if(a != NULL){
      freeNode(a->next);
      free(a);
    }
}
//free a graph
void freeGraph(graph *g){
  for(int i = 0; i < g->nr_size; i++){
    if(!isEmpty((*g->gr)[(*g->nodes)[i]]))
      freeNode((*g->gr)[(*g->nodes)[i]].next);
  }
  free(g->gr);
  free(g->nodes);
  free(g);
}
//depth first algorithm for calculating the depth of the network
int dfs(graph *g, node *n, bool vi[]){

  if(n->value == 0)
  {
    vi[n->value] = true;
    return 1;
  }
  if(vi[n->value] == true)
    return -1;
  vi[n->value] = true;
  node *a = n;
  int x = 0;
  while(a != NULL && x != -1){

    int i = dfs(g, &(*g->gr)[a->value], vi);
    if(a->value > 0)
      vi[a->value] = true;
    if(i >= 0){
      if(x < i)
        x = i;
    }
    else
      x = -1;
    a = a -> next;
  }
  if(x != -1)
    return 1 + x;
  return -1;
}

int calculateTreeDepth(graph *g){
  unsigned int n = g->r_size + 1;
  bool vi[g->size + 2];
  for(int i = 0; i < n; i++){
    vi[i] = false;
  }
  vi[(*g->nodes)[0]] = true;
  int rez = dfs(g, &(*g->gr)[(*g->nodes)[0]], vi);
  for(int i = 0; i < g->nr_size ; i++){
    if((*g->nodes)[i] != -1)
      if(vi[(*g->nodes)[i]] == false)
        return -1;
      }
  return rez;
}

//print on screen the given graph
void show(graph *g){
  for(int i = 0; i < g->nr_size;i++){
      printf("%d ->",(*g->nodes)[i]);
      node *a = &(*g->gr)[(*g->nodes)[i]];
      if(!isEmpty(*a)){
        while(a != NULL){
          printf(" %d distance = %lf  ", a->value, a->w);
          a = a->next;
        }
      }
      printf("\n");
  }
}
bool chk(int rez[], int i){
  for(int j = 0; j < i - 2; j++)
    for(int k = j + 1; k < i-1; k++)
      if(rez[j] == rez[k])
        return false;
  if(rez[0] == rez[i-1])
    return true;
  return false;
}

//depth first algorithm for calculating the hamiltonian cycles
void dfs2(graph *g, node *n, int rez[], int i, bool *ok, double dist){
  if(i == g->nr_size + 1){
    if(chk(rez, i) == true){
      for(int j = 0; j < i; j++){
        printf("%d", rez[j]);
        if(j != i - 1)
          printf("->");
      }
      printf ("   distance = %lf\n",dist);
      *ok = true;
    }
  }
  else if(i < g->nr_size + 1){
    node *a = n;
    while(a != NULL){
      rez[i] = a->value;
      dist += a->w;
      dfs2(g, &(*g->gr)[a->value], rez, i + 1, ok, dist);
      dist -= a->w;
      a = a->next;
    }
  }
}
void hamilt(graph *g){
  int rez[g->nr_size + 1];
  bool ok = 0;
  rez[0] = 1;
  dfs2(g, &(*g->gr)[1], rez, 1, &ok, 0);
  if(ok == false)
    printf("There are no hamiltonian cycles in the given graph\n");
}

//creates a matrix from the given network
void createMatrix(graph *g, double m[][g->r_size + 1]){
  for(int i = 0; i < g->r_size + 1; i++)
    for(int j = 0; j < g->r_size + 1; j++)
      m[i][j] = 0;
  for(int i = 0; i < g->nr_size; i++){
    node* n = &(*g->gr)[(*g->nodes)[i]];
    while(n != NULL){
      m[(*g->nodes)[i]][n->value] = n->w;
      n = n->next;
    }
  }
}
// roy Roy-Floyd algorithm on a given matrix
void royFloyd(int n, double m[][n]){
  for(int k = 1; k < n; k++)
    for(int i = 1; i < n; i++)
      for(int j = 1; j < n; j++)
        if(m[i][k] && m[k][j] && (m[i][j] > m[i][k] + m[k][j] || !m[i][j]) && i != j)
          m[i][j] = m[i][k] + m[k][j];
}


void showMatrix(int n, double m[][n]){
  for(int i = 1; i < n; i++){
      for(int j = 1; j < n; j++)
        printf("%10lf ",m[i][j]);
      printf("\n");
    }
}
//depth first algorithm for strongly connected algorithm
void dfs3(graph *g, node *n, bool vi[]){
  if(vi[n->value] != true){
    vi[n->value] = true;
    node *a = &(*g->gr)[n->value];
    while(a != NULL){
      dfs3(g, a, vi);
      a = a->next;
    }
  }
}
//check if a graph is strongly connected
bool isStronglyConnected(graph *g){
  bool vi[g->size];
  for(int j = 0; j < g->nr_size; j++){
    for(int i = 0 ;i < g->size; i++)
      vi[i] = false;
    dfs3(g, &(*g->gr)[(*g->nodes)[j]], vi);
    for(int i = 0; i < g->nr_size; i++)
      if(vi[(*g->nodes)[i]] == false)
        return false;
  }
  return true;
}
//check if a graph is connected
bool isConnected(graph *g){
  if(g->nr_size != 1){
    for(int i = 0; i < g->nr_size; i++){
      node *a = &(*g->gr)[(*g->nodes)[i]];
      while(a != NULL){
        push(g, a->value, (*g->nodes)[i], -1);
        a = a->next;
      }
    }
  }
  return isStronglyConnected(g);
}
//search for an eulerian cycle
void euler(graph *g, int a, int v[], long n, long *mx){
  if(n > *mx)
    *mx = n;
  node *m = &(*g->gr)[a];
  while(m!= NULL){
    if(m->w != -1){
      m->w = -1;
      euler(g, m->value, v, n + 1, mx);
    }
    m = m->next;
  }
  v[n] = a;
}

bool isDegree(graph *g){
  pair a[g->r_size + 1];
  for(int i = 0 ; i < g->nr_size ; i++)
    a[(*g->nodes)[i]] = (pair){0,0};
  for(int i = 0 ; i < g->nr_size ; i++){
    node *n = & (*g->gr)[(*g->nodes)[i]];
    while (n!= NULL){
      a[(*g->nodes)[i]].out ++;
      a[n->value].in ++;
      n = n->next;
    }
  }
for(int i = 0 ; i < g->nr_size ; i++){
    if(a[(*g->nodes)[i]].in != a[(*g->nodes)[i]].out)
      return false;
}
  return true;
}

void isEuler(graph *g){
  int v[1l * g->r_size * g->r_size + 1];
  long n = 0;
  if(isStronglyConnected(g) == 1 && isDegree(g) == 1){
      euler(g, (*g->nodes)[0], v, 0, &n);
      for(int i = 0; i <= n; i++){
        if(i != 0)
          printf("->");
        printf("%d",v[i]);
      }
      printf("\n");
    }
  else{
    printf("The given graph does not have an eulerian cycle\n");
  }
}

//algorithm for the maximal clique in the graph
void BronKerbosch(graph* g, int r[], int p[], int x[], int n, int m, int b, int rez[], int *srez){
  if(b == 0 && m == 0){
    if(*srez < n){
      for(int i = 0; i < n; i++)
        rez[i] = r[i];
        *srez = n;
      }
  }

  for(int i = 0; i < m; i++){
    if(p[i] != -1){
      r[n] = p[i];
      n++;
      int pn[2 * g->nr_size + 1];
      int j = 0;
      int xn[2 * g->nr_size + 1];
      int k = 0;
      for(int d = 0; d < m; d++){
        node *n = &(*g->gr)[p[i]];
        while(n != NULL){
          if(n->w != -1){
            if(n->value == p[d]){
              pn[j] = p[d];
              j++;
            }
          }
          n = n->next;
        }
      }
      for(int d = 0; d < b; d++){
        node *n = &(*g->gr)[p[i]];
        while(n != NULL){
          if(n->w != -1){
            if(n->value == x[d]){
              xn[k] = x[d];
              k++;
            }
          }
          n = n->next;
        }
      }
      BronKerbosch(g, r, pn, xn, n, j, k, rez, srez);
      n--;
      p[i] = -1;
      bool ok = 0;
      for(int d = 0; d < b; d++)
        if(x[d] == p[i])
          ok = 1;
        if(ok == 0){
        x[b] = p[i];
        b++;
      }
      i--;
    }
  }
}

bool isCon(graph *g, int x, int i){
  node *a = &(*g->gr)[x];
  while(a != NULL){
    if(a->value == i)
      return true;
    a = a->next;
  }
  return false;
}

void isBronKerbosch(graph *g){
  for(int i = 0; i < g->nr_size; i++){
    node *a = &(*g->gr)[(*g->nodes)[i]];
    while(a != NULL){
      if(!isCon(g, a->value, (*g->nodes)[i])){
        a->w = -1;
      }
      a = a->next;
    }
  }
  int rez[2 * g->nr_size + 1];
  int r[2 * g->nr_size + 1], p[2 * g->nr_size + 1], x[2 * g->nr_size + 1];
  int n = 0, m = 0, b = 0, srez = 0;
  for(int i = 0; i < g->nr_size; i++)
    p[i] = (*g->nodes)[i];
  m = g ->nr_size;
  BronKerbosch(g, r, p, x, n, m, b, rez, &srez);
  if(srez != 1){
    printf("Maximal clique is : ");
    for(int i = 0; i < srez; i++)
      printf("%d ",rez[i]);
    printf("\n");
  }
  else printf("There are no cliques in the graph\n");
}

//chose which function to use
void fnx(const char s[], graph *g){
  if(strcmp(s, "hamiltonian") == 0){
    if(g->r_size > 9000000)
      printf("The nodes are to big to calculate the eulerian cycle\n");
    else
      hamilt(g);
  }
  else if(strcmp(s, "show") == 0)
    show(g);
  else if(strcmp(s, "matrix") == 0){
    if(g->r_size > 1000)
      printf("The values of the nodes are too big for the matrix to be calculated\n");
    else{
      double m[g->r_size + 1][g->r_size + 1];
      createMatrix(g, m);
      if(g->r_size > 100){
        printf("the matrix was created but it will not be printed because of its size\n");
      }
      else
        showMatrix(g->r_size + 1, m);
      }
  }
  else if(strcmp(s,"distance") == 0){
    if(g->r_size > 700)
      printf("The values of the nodes are too big for the distance matrix to be calculated\n");
    else{
      double m[g->r_size + 1][g->r_size + 1];
      createMatrix(g, m);
      royFloyd(g->r_size + 1, m);
      if(g->r_size > 100){
        printf("the matrix with the distances was created but it will not be printed because of its size\n");
      }
      else
        showMatrix(g->r_size + 1, m);
    }
  }
  else if(strcmp(s,"sconnected") == 0){
    if(isStronglyConnected(g))
      printf("The graph is strongly connected\n");
    else
      printf("The graph is not strongly connected\n");
  }
  else if(strcmp(s,"connected") == 0){
    if(isConnected(g))
      printf("The graph is connected\n");
    else
      printf("The graph is not connected\n");
  }
  else if(strcmp(s,"eulerian") == 0){
    if(g->r_size > 1000)
      printf("The nodes are to big to calculate the eulerian cycle\n");
    else
      isEuler(g);
  }
  else if (strcmp(s,"clique") == 0){
    if(g->r_size > 9000000)
      printf("The nodes are to big to calculate the cliques\n");
    else
      isBronKerbosch(g);
  }
  else printf("Invalid function\n");
}

bool verif(graph * g, char s[]){
  int n = strlen(s);
  for(int i = 0; i < n; i++)
    if(s[i] < '0' || s[i] > '9')
      return false;
  if(!isNode(g,atoi(s)))
    return false;
  return true;
}

void dijkstra(int i, graph *g, double dis[]){
  bool ap[g->r_size + 1];
  for(int i = 0; i < g->r_size + 1; i++){
    ap[i] = false;
    dis[i] = DBL_MAX;
  }
  dis[i] = 0;
  int hp[g->nr_size + 1];
  int n = 0;
  hp[n] = i;
  n++;
  while(n != 0){
    int nod = hp[n - 1];
    ap[nod] = 0;
    n--;
    node *a = &(*g->gr)[nod];
    while(a != NULL){
      if(dis[a->value] > dis[nod] + a->w){
        dis[a->value] = dis[nod] + a->w;
        if(ap[a->value] == false){
          hp[n] = a->value;
          n++;
          ap[a->value] = 1;
        }
      }
      a = a->next;
    }
  }
}

void isDijs(char a1[], char a2[], graph *g){
  if(!verif(g, a2)){
    printf("Invalid Node!\n");
    exit(1);
  }
  if(strcmp(a1, "distance") == 0){
    double dis[g->r_size + 1];
    dijkstra(atoi(a2), g, dis);
    for(int i = 0; i < g->nr_size; i++)
      printf("   %d       ", (*g->nodes)[i]);
    printf("\n");
    for(int i = 0; i < g->nr_size; i++){
      if(dis[(*g->nodes)[i]] == DBL_MAX)
        printf("  n/a    ");
      else
        printf("%10lf ", dis[(*g->nodes)[i]]);
    }
    printf("\n");
  }
  else
    printf("Invalid function \n");
}

void assert(int line, bool b) {
    if (b) return;
    printf("The test on line %d fails.\n", line);
    exit(1);
}

bool eq(int a[], int b[], int n){
  for(int i = 0; i < n ; i++)
    if(a[i] != b[i])
      return false;
  return true;
}

bool gf(graph *g ,int a[], int n){
  int x[g->size];
  int m = 0;
  for(int i = 0; i < n; i++){
    int ok = 0;
    for(int i = 0; i < m; i++)
      if(x[i] == (*g->nodes)[i])
        ok = 1;
    if(ok == 0){
      x[m] = (*g->nodes)[i];
      m++;
    }
    node *n = &(*g->gr)[(*g->nodes)[i]];
    while(n != NULL){
      if(n->value != 0){
        int ok = 0;
        for(int i = 0; i < m; i++)
          if(x[i] == n->value)
            ok = 1;
        if(ok == 0){
          x[m] = n->value;
          m++;
        }
      }
      n = n->next;
    }
  }
  qsort(x, m, sizeof(int), comp);
  return eq(x, (*g->nodes), m);
}
void testCreateNetwork(){
  char s[] = "1-2/3,1-3/2.4,1-4";
  graph *g = createNetwork(s);
  assert(__LINE__, g->size == 24);
  assert(__LINE__, g->r_size == 4);
  assert(__LINE__, g->nsize == 24);
  assert(__LINE__, g->nr_size == 4);
  int te[] = {1,2,3,4};
  assert(__LINE__, (*g->gr)[1].w == 3);
  assert(__LINE__, (*g->gr)[1].next->w == 2.4);
  assert(__LINE__, (*g->gr)[1].next->next->w == 1);
  assert(__LINE__, eq(te, (*g->nodes), g->nr_size));
  assert(__LINE__, gf(g, (*g->nodes), g->nr_size));
  strcpy(s, "1-1000");
  freeGraph(g);
  g = createNetwork(s);
  assert(__LINE__, g->size == 1369);
  assert(__LINE__, g->r_size == 1000);
  assert(__LINE__, g->nsize == 24);
  assert(__LINE__, g->nr_size == 2);
  int te2[] ={1,1000};
  assert(__LINE__, eq(te2, (*g->nodes), g->nr_size));
  assert(__LINE__, gf(g, (*g->nodes), g->nr_size));
  freeGraph(g);
}

void testDepth(){
  char s[100];
  strcpy(s,"1-2/3,1-3/2.4,1-4");
  graph *g = createNetwork(s);
  assert(__LINE__, calculateTreeDepth(g) == 2);
  freeGraph(g);
  strcpy(s, "1-2,2-3,3-1");
  g = createNetwork(s);
  assert(__LINE__, calculateTreeDepth(g) == -1);
  freeGraph(g);
  strcpy(s, "1-2,1-3,1-6,6-7,7-4");
  g = createNetwork(s);
  assert(__LINE__, calculateTreeDepth(g) == 4);
  freeGraph(g);
  strcpy(s, "1-2,3");
  g = createNetwork(s);
  assert(__LINE__, calculateTreeDepth(g) == -1);
  freeGraph(g);
  strcpy(s, "1-1,1-2");
  g = createNetwork(s);
  assert(__LINE__, calculateTreeDepth(g) == -1);
  freeGraph(g);
  strcpy(s, "1-1000,1000-2");
  g = createNetwork(s);
  assert(__LINE__, calculateTreeDepth(g) == 3);
  freeGraph(g);
  strcpy(s, "1-10000,10000-200,2000-3");
  g = createNetwork(s);
  assert(__LINE__, calculateTreeDepth(g) == -1);
  freeGraph(g);
}

bool eqM(int n, double a[][n], double b[][n]){
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      if(a[i][j] != b[i][j])
        return false;
  return true;
}

void testCreateMatrix(){
  double a[6][6] = {{0,0,0,0,0,0},{0,0,3,9,8,3},{0,5,0,1,4,2},{0,6,6,0,4,5},{0,2,9,2,0,7},{0,7,9,3,2,0}};
  double c[4][4] = {{0,0,0,0},{0,0,1,1},{0,0,0,0},{0,0,0,0}};
  double b[100][100];
  char s[1000];
  strcpy(s,"1-2,1-3");
  graph *g = createNetwork(s);
  createMatrix(g,b);
  assert(__LINE__, eqM(4, c, b));
  freeGraph(g);
  strcpy(s, "1-2/3,1-3/9,1-4/8,1-5/3,2-1/5,2-3/1,2-4/4,2-5/2,3-1/6,3-2/6,3-4/4,3-5/5,4-1/2,4-2/9,4-3/2,4-5/7,5-1/7,5-2/9,5-3/3,5-4/2");
  g = createNetwork(s);
  createMatrix(g,b);
  assert(__LINE__, eqM(6, a, b));
  freeGraph(g);
}

void testRoyFloyd(){
  double a[6][6] = {{0,0,0,0,0,0},{0,0,3,4,5,3},{0,5,0,1,4,2},{0,6,6,0,4,5},{0,2,5,2,0,5},{0,4,7,3,2,0}};
  char s[1000];
  double b[100][100];
  strcpy(s,"1-2,1-3,3-1");
  double c[4][4] = {{0,0,0,0},{0,0,1,1},{0,0,0,0},{0,1,2,0}};
  graph *g = createNetwork(s);
  createMatrix(g, b);
  royFloyd(g->r_size + 1, b);
  assert(__LINE__, eqM(4, b, c));
  freeGraph(g);
  strcpy(s, "1-2/3,1-3/9,1-4/8,1-5/3,2-1/5,2-3/1,2-4/4,2-5/2,3-1/6,3-2/6,3-4/4,3-5/5,4-1/2,4-2/9,4-3/2,4-5/7,5-1/7,5-2/9,5-3/3,5-4/2");
  g = createNetwork(s);
  createMatrix(g, b);
  royFloyd(g->r_size + 1, b);
  assert(__LINE__, eqM(6, a, b));
  freeGraph(g);
}

void testSconnected(){
  char s[1000];
  strcpy(s, "1-2,2-300,300-1");
  graph *g = createNetwork(s);
  assert(__LINE__, isStronglyConnected(g) == true);
  freeGraph(g);
  strcpy(s, "2");
  g = createNetwork(s);
  assert(__LINE__, isStronglyConnected(g) == false);
  freeGraph(g);
  strcpy(s, "2000-2000");
  g = createNetwork(s);
  assert(__LINE__, isStronglyConnected(g) == true);
  freeGraph(g);
  strcpy(s, "2-2,3");
  g = createNetwork(s);
  assert(__LINE__, isStronglyConnected(g) == false);
  freeGraph(g);
  strcpy(s, "1-2,3-4,2-3");
  g = createNetwork(s);
  assert(__LINE__, isStronglyConnected(g) == false);
  freeGraph(g);
  strcpy(s, "4-5,5-4,5-800,800-7,7-6,6-7,6-5");
  g = createNetwork(s);
  assert(__LINE__, isStronglyConnected(g) == true);
  freeGraph(g);
}

void testConnected(){
  char s[1000];
  strcpy(s, "1-2,2-3,3-1");
  graph *g = createNetwork(s);
  assert(__LINE__, isConnected(g) == true);
  freeGraph(g);
  strcpy(s, "2");
  g = createNetwork(s);
  assert(__LINE__, isConnected(g) == false);
  freeGraph(g);
  strcpy(s, "200-200");
  g = createNetwork(s);
  assert(__LINE__, isConnected(g) == true);
  freeGraph(g);
  strcpy(s, "1-2");
  g = createNetwork(s);
  assert(__LINE__, isConnected(g) == true);
  freeGraph(g);
  strcpy(s, "4-5,5-4,5-1000,1000-7,7-6,6-7,6-5");
  g = createNetwork(s);
  assert(__LINE__, isConnected(g) == true);
  freeGraph(g);
  strcpy(s, "4-5,5-8,8-7,7-6,6-5");
  g = createNetwork(s);
  assert(__LINE__, isConnected(g) == true);
  freeGraph(g);
  strcpy(s, "1-2,1-4,2-3,4-8,8-7,7-6,6-5");
  g = createNetwork(s);
  assert(__LINE__, isConnected(g) == true);
  freeGraph(g);
}

void testEulerian(){
  char s[1000];
  int v[1000];
  long n = 0;
  graph *g;
  strcpy(s, "1-2,2-3,3-4,4-1");
  g = createNetwork(s);
  euler(g, (*g->nodes)[0], v, 0, &n);
  int rez1[] = {1,2,3,4,1};
  assert(__LINE__, isStronglyConnected(g) == 1 && isDegree(g) == 1);
  assert(__LINE__, eq(v ,rez1 ,n));
  freeGraph(g);
  strcpy(s, "2-2");
  int rez2[] = {2,2};
  g = createNetwork(s);
  euler(g, (*g->nodes)[0], v, 0, &n);
  assert(__LINE__, isStronglyConnected(g) == 1 && isDegree(g) == 1);
  assert(__LINE__, eq(v ,rez2 ,2));
  freeGraph(g);
  strcpy(s, "20-200,200-20");
  int rez3[] = {20,200,20};
  g = createNetwork(s);
  euler(g, (*g->nodes)[0], v, 0, &n);
  assert(__LINE__, isStronglyConnected(g) == 1 && isDegree(g) == 1);
  assert(__LINE__, eq(v ,rez3 ,3));
  freeGraph(g);
  strcpy(s, "1-2,1-3,2-2,2-3,3-4,3-4");
  g = createNetwork(s);
  euler(g, (*g->nodes)[0], v, 0, &n);
  assert(__LINE__, isStronglyConnected(g) == 0 && isDegree(g) == 0);
  freeGraph(g);
  strcpy(s, "1-2,1-3,2-2,2-3,3-4,2-1,3-1,3-2,4-3");
  int rez[] = {1,2,2,3,4,3,1,3,2,1};
  g = createNetwork(s);
  assert(__LINE__, isStronglyConnected(g) == 1 && isDegree(g) == 1);
  euler(g, (*g->nodes)[0], v, 0, &n);
  assert(__LINE__, eq(v ,rez ,n));
  freeGraph(g);
}

void testClique(){
  graph *g;
  char s[1000];
  strcpy(s, "1-2,1-3,2-3");
  g = createNetwork(s);
  for(int i = 0; i < g->nr_size; i++){
    node *a = &(*g->gr)[(*g->nodes)[i]];
    while(a != NULL){
      if(!isCon(g, a->value, (*g->nodes)[i])){
        a->w = -1;
      }
      a = a->next;
    }
  }
  int rez[2 * g->nr_size + 1];
  int r[2 * g->nr_size + 1], p[2 * g->nr_size + 1], x[2 * g->nr_size + 1];
  int n = 0, m = 0, b = 0, srez = 0;
  BronKerbosch(g, r, p, x, n, m, b, rez, &srez);
  assert(__LINE__, srez == 0);
  freeGraph(g);
  strcpy(s, "1-2,2-1,1-5,5-1,2-5,5-2,2-3,3-2,3-4,4-3,4-5,5-4,4-6,6-4");
  g = createNetwork(s);
  for(int i = 0; i < g->nr_size; i++){
    node *a = &(*g->gr)[(*g->nodes)[i]];
    while(a != NULL){
      if(!isCon(g, a->value, (*g->nodes)[i])){
        a->w = -1;
      }
      a = a->next;
    }
  }
  n = 0, m = 0, b = 0, srez = 0;
  BronKerbosch(g, r, p, x, n, m, b, rez, &srez);
  int rez1[] = {1,2,5};
  assert(__LINE__, eq(rez, rez1, srez));
  freeGraph(g);
}

bool eqD(double a[], double b[], int n){
  for(int i = 0; i < n; i++)
    if(a[i] != b[i])
      return false;
  return true;
}
void testDij(){
  graph *g;
  char s[1000];
  strcpy(s, "1-2/1,1-4/2,4-3/4,2-3/2,4-5/3,3-5/6");
  g = createNetwork(s);
  double rez[] = {5,0,1,3,2};
  double dis[1000];
  dijkstra(1, g, dis);
  assert(__LINE__, eqD(rez,dis,5));
  dijkstra(2, g, dis);
  double rez2[] = {8,DBL_MAX,0,2,DBL_MAX};
  assert(__LINE__,eqD(rez2,dis,5));
  freeGraph(g);
}

void test(){
  testCreateNetwork();
  testDepth();
  testCreateMatrix();
  testRoyFloyd();
  testSconnected();
  testConnected();
  testEulerian();
  testClique();
  testDij();
  printf("All tests run OK!\n");
}

int main(int n, char *args[n]) {
    setbuf(stdout, NULL);
    if (n == 1) {
        test();
    }
    else if (n == 2) {
      graph *g = createNetwork(args[1]);
      printf("number of nodes: %d\n", g->nr_size);
      int x = calculateTreeDepth(g);
      if(x < 0)
        printf("the given graph is not a tree\n");
      else
        printf("tree depth: %d\n", x);
      freeGraph(g);

    }
    else if(n==3){
      graph *g = createNetwork(args[1]);
      fnx(args[2], g);
      freeGraph(g);
    }
    else if(n==4){
      graph *g = createNetwork(args[1]);;
      isDijs(args[2],args[3],g);
      freeGraph(g);
    }
    else {
        fprintf(stderr, "Invalid number of arguments\n");
        exit(1);
    }
    return 0;
}
