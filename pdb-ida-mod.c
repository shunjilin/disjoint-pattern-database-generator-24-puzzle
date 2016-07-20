/* This program performs iterative-deepening A* on the 24 puzzles, using
  disjoint pattern database heuristic functions. This version solves the twenty
  puzzle with four pattern databases, splitting the tiles into groups of 6.
  This version uses the pattern table containing all the moves, and doesn't
  compute manhattan distances.  This version uses direct-access pattern
  databases, without the n^2 hashing function.  This version computes the
  maximum of the original pattern databases and their reflection about the main
  diagonal.  This version uses only one copy of the 2x3 pattern, rotating the
  puzzle to get the values for the other 2 copies. This version only computes
  the reflected heuristic function if the original function plus the depth
  doesn't exceed the threshold. This is normally the fastest version. Written
  by Richard E. Korf, Computer Science Department, University of California,
  Los Angeles, CA 90095. */

#include <stdio.h>                                    /* standard I/O library */

#define NUMBER 100                /* number of problem instances to be solved */
#define X 5                                     /* squares in the x dimension */
#define SIZE 25                                    /* total number of squares */
#define MAXMOVES 125                               /* maximum number of moves */
#define TABLESIZE 244140625   /* bytes in direct-access database array (25^6) */

FILE *infile;                              /* pointer to heuristic table file */

int s[SIZE];     /* state of puzzle: for each position, tile in that position */
int inv[SIZE];    /* inverse state: for each tile, its position in the puzzle */

struct operators
{int num;                                  /* number of applicable oprs: 2..4 */
 int pos[4];} oprs[SIZE];     /* position of adjacent tiles for each position */

unsigned char h0[TABLESIZE];        /* heuristic tables for pattern databases */
unsigned char h1[TABLESIZE];

/* tiles in each regular pattern */
/* {1,2,5,6,7,12} {3,4,8,9,13,14} {10,11,15,16,20,21} {17,18,19,22,23,24} */

/* tiles in each reflected pattern, in the same order as above*/
/* {5,10,1,6,11,12} {15,20,16,21,17,22} {2,7,3,8,4,9} {13,18,23,14,19,24} */


/* the pattern that each tile is in */
int whichpat[25] = {0,0,0,1,1,0,0,0,1,1,2,2,0,1,1,2,2,3,3,3,2,2,3,3,3};

/* the reflected pattern that each tile is in */
int whichrefpat[25] = {0,0,2,2,2,0,0,2,2,2,0,0,0,3,3,1,1,1,3,3,1,1,1,3,3};

/* the position of each tile in order, reflected about the main diagonal */
int ref[SIZE] = {0,5,10,15,20,1,6,11,16,21,2,7,12,17,22,3,8,13,18,23,4,9,14,19,24};

/* rotates the puzzle 90 degrees */
int rot90[25] = {20,15,10,5,0,21,16,11,6,1,22,17,12,7,2,23,18,13,8,3,24,19,14,9,4};

/* composes the reflection and 90 degree rotation into a single array */
int rot90ref[25] = {20,21,22,23,24,15,16,17,18,19,10,11,12,13,14,5,6,7,8,9,0,1,2,3,4};

/* rotates the puzzle 180 degrees */
int rot180[25] = {24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0};

/* composes the reflection and 180 degree rotation into a single array */
int rot180ref[25] = {24,19,14,9,4,23,18,13,8,3,22,17,12,7,2,21,16,11,6,1,20,15,10,5,0};

int path[MAXMOVES];                                          /* solution path */

int thresh;                                        /* search cutoff threshold */
double generated;                 /* number of states generated per iteration */
double total;                             /* total number of states generated */
double grandtotal;                 /* total nodes generated in all iteration s*/

/* INITOPS initializes the operator table. */

initops ()

{int blank;                                   /* possible positions of blank */

 for (blank = 0; blank < SIZE; blank++)  /* for each possible blank position */
  {oprs[blank].num = 0;                               /* no moves initially */
   if (blank > X - 1)                                       /* not top edge */
     oprs[blank].pos[oprs[blank].num++] = blank - X;       /* add a move up */
   if (blank % X > 0)                                      /* not left edge */
     oprs[blank].pos[oprs[blank].num++] = blank - 1;     /* add a move left */
   if (blank % X < X - 1)                                 /* not right edge */
     oprs[blank].pos[oprs[blank].num++] = blank + 1;    /* add a move right */
   if (blank < SIZE - X)                                 /* not bottom edge */
     oprs[blank].pos[oprs[blank].num++] = blank + X;}}   /* add a move down */

/* INPUT accepts an initial state from the terminal, assuming it is preceded by
a problem number. It stores it in the state vector and returns the position of
the blank tile. */

input (s)

int s[SIZE];                                                  /* state vector */

{int index;                                        /* index to tile positions */
 int blank;                                         /* position of blank tile */

 //scanf ("%*d");                                   /* skip over problem number */
 for (index = 0; index < SIZE; index++)                  /* for each position */
  {scanf ("%d", &s[index]);                   /* input tile in that position */
   if (s[index] == 0) blank = index;}      /* note blank position in passing */
 return (blank);}

/* HASH0 takes an inverse state, and maps the tiles in the 0 pattern to an
  integer that represents those tile positions uniquely.  It then returns the
  actual heuristic value from the pattern database. */

unsigned int hash0 ()

{int hashval;                                   /* index into heuristic table */

 hashval = ((((inv[1]*SIZE+inv[2])*SIZE+inv[5])*SIZE+inv[6])*SIZE+inv[7])*SIZE+inv[12];
 return (h0[hashval]);}                       /* total moves for this pattern */

/* HASHREF0 takes an inverse state, and maps the tiles in the 0 pattern to an
  integer that represents those tile positions uniquely.  It then returns the
  actual heuristic value from the pattern database.  This is the reflected
  version of this database. */


unsigned int hashref0 ()

{int hashval;                                   /* index into heuristic table */

 hashval = (((((ref[inv[5]] * SIZE
               + ref[inv[10]]) * SIZE
              + ref[inv[1]]) * SIZE
             + ref[inv[6]]) * SIZE
            + ref[inv[11]]) * SIZE
           + ref[inv[12]]);

 return (h0[hashval]);}                       /* total moves for this pattern */

/* HASH1 takes an inverse state, and maps the tiles in the 1 pattern to an
  integer that represents those tile positions uniquely.  It then returns the
  actual heuristic value from the pattern database. */

unsigned int hash1 ()

{int hashval;                                   /* index into heuristic table */

 hashval = ((((inv[3]*SIZE+inv[4])*SIZE+inv[8])*SIZE+inv[9])*SIZE+inv[13])*SIZE+inv[14];

 return (h1[hashval]);}                       /* total moves for this pattern */

/* HASHREF1 takes an inverse state, and maps the tiles in the 1 pattern to an
  integer that represents those tile positions uniquely.  It then returns the
  actual heuristic value from the pattern database.  This is the reflected
  version of this database. */

unsigned int hashref1 ()

{int hashval;                                   /* index into heuristic table */

 hashval = (((((ref[inv[15]] * SIZE
               + ref[inv[20]]) * SIZE
              + ref[inv[16]]) * SIZE
             + ref[inv[21]]) * SIZE
            + ref[inv[17]]) * SIZE
           + ref[inv[22]]);

 return (h1[hashval]);}                       /* total moves for this pattern */

/* HASH2 takes an inverse state, and maps the tiles in the 2 pattern to an
  integer that represents those tile positions uniquely.  It then returns the
  actual heuristic value from the pattern database. */

unsigned int hash2 ()

{int hashval;                                   /* index into heuristic table */

 hashval = ((((rot180[inv[21]] * SIZE
              + rot180[inv[20]]) * SIZE
             + rot180[inv[16]]) * SIZE
            + rot180[inv[15]]) * SIZE
           + rot180[inv[11]]) * SIZE
          + rot180[inv[10]];

 return (h1[hashval]);}                       /* total moves for this pattern */

/* HASHREF2 takes an inverse state, and maps the tiles in the 2 pattern to an
  integer that represents those tile positions uniquely.  It then returns the
  actual heuristic value from the pattern database.  This is the reflected
  version of this database. */

unsigned int hashref2 ()

{int hashval;                                   /* index into heuristic table */

 hashval = (((((rot180ref[inv[9]] * SIZE
               + rot180ref[inv[4]]) * SIZE
              + rot180ref[inv[8]]) * SIZE
             + rot180ref[inv[3]]) * SIZE
            + rot180ref[inv[7]]) * SIZE
           + rot180ref[inv[2]]);

 return (h1[hashval]);}                       /* total moves for this pattern */

/* HASH3 takes an inverse state, and maps the tiles in the 3 pattern to an
  integer that represents those tile positions uniquely.  It then returns the
  actual heuristic value from the pattern database. This version rotates the
  puzzle and uses the 1 pattern database. */

unsigned int hash3 ()

{int hashval;                                   /* index into heuristic table */

 hashval = ((((rot90[inv[19]] * SIZE
              + rot90[inv[24]]) * SIZE
             + rot90[inv[18]]) * SIZE
            + rot90[inv[23]]) * SIZE
           + rot90[inv[17]]) * SIZE
          + rot90[inv[22]];

 return (h1[hashval]);}                       /* total moves for this pattern */

/* HASHREF3 takes an inverse state, and maps the tiles in the 3 pattern to an
  integer that represents those tile positions uniquely.  It then returns the
  actual heuristic value from the pattern database.  This is the reflected
  version of this database. */

unsigned int hashref3 ()

{int hashval;                                   /* index into heuristic table */

 hashval = (((((rot90ref[inv[23]] * SIZE
               + rot90ref[inv[24]]) * SIZE
              + rot90ref[inv[18]]) * SIZE
             + rot90ref[inv[19]]) * SIZE
            + rot90ref[inv[13]]) * SIZE
           + rot90ref[inv[14]]);

 return (h1[hashval]);}                       /* total moves for this pattern */

/* SEARCH performs one depth-first iteration of the search, cutting off
 when the depth plus the heuristic evaluation exceeds THRESH. If it succeeds,
 it returns 1 and records the sequence of tiles moved in the solution.
 Otherwise, it returns 0 */

search (blank, oldblank, g, add0, add1, add2, add3, addr0, addr1, addr2, addr3)

int blank;                                      /* current position of blank */
int oldblank;                                  /* previous position of blank */
int g;                                            /* current depth of search */
int add0, add1, add2, add3;   /* additional moves from each regular database */
int addr0, addr1, addr2, addr3; /* additional moves from each reflected database */

{int index;                                     /* index into operator array */
 int newblank;                                /* blank position in new state */
 int tile;                                               /* tile being moved */
 int nadd0, nadd1, nadd2, nadd3; /*new additional moves from regular database*/
 int naddr0, naddr1, naddr2, naddr3; /*new additional moves from reflected database*/
 int nadd;               /* total new additional moves from regular database */
 int naddr;            /* total new additional moves from reflected database */

 for (index = 0; index < oprs[blank].num; index++)    /* each applicable op */
  if ((newblank = oprs[blank].pos[index]) != oldblank) /*not inv last move */
    {tile = s[newblank];            /* tile moved is in new blank position */
    s[blank] = tile;                                   /* make actual move */
    inv[tile] = blank;                   /* maintain inverse state as well */
    generated++;                                  /* count nodes generated */

    if (whichpat[tile] == 0)                       /* tile is in 0 pattern */
      {nadd0 = hash0();                 /* additional moves from 0 pattern */
	/*if ((nadd0 + 1) - add0 < 0 ) {
	  printf("diff is %d\n", nadd0 + 1 - add0);
	  printf("tile moved is %d\n", tile);
	  printf("nadd0 is %d\n", nadd0);
	  printf("add0 is %d\n", add0);
	  printf("old blank is %d\n", oldblank);
	  printf("current blank is %d\n", blank);
	  printf("new blank is %d\n", newblank);
	  printf("-------------------------------");
	  printf("\n");
	  for (int i = 0; i < 5; ++ i) {
	    printf("%3d", s[i]);
	  }
	  printf("\n");
	  for (int i = 5; i < 10; ++i) {
	    printf("%3d", s[i]);
	  }
	  printf("\n");
	  for (int i = 10; i < 15; ++i) {
	    printf("%3d", s[i]);
	  }
	  printf("\n");
	  for (int i = 15; i < 20; ++i) {
	    printf("%3d", s[i]);
	  }
	  printf("\n");
	  for (int i = 20; i < 25; ++i) {
	    printf("%3d", s[i]);
	  }
	  printf("\n");
	  printf("nadd0 is %d\n", nadd0);
	  printf("add1 is %d\n", add1);
	  printf("add2 is %d\n", add2);
	  printf("add3 is %d\n", add3);
	  }*/
      nadd = nadd0 + add1 + add2 + add3;     /* new moves from reg pattern */
      if (nadd + g >= thresh) goto next; /*already > thresh, skip this move*/

      if (whichrefpat[tile] == 0)        /* tile is in 0 reflected pattern */
        {naddr0 = hashref0();           /* additional moves from 0 pattern */
	  //if ((naddr0 + 1) - addr0 < 0 ) {
	      //printf("addr0 inconsistent\n");
	    //}
        naddr = naddr0 + addr1 + addr2 + addr3;  /* total additional moves */
        if (naddr + g < thresh)                 /* naddr + g + 1 <= thresh */
          if ((nadd == 0) ||                   /* goal state is reached or */
              (search(newblank,blank,g+1,nadd0,add1,add2,add3,naddr0,addr1,addr2,addr3)))
            {path[g] = tile;                          /* record tile moved */
            return (1);}}                             /* exit with success */
      else                               /* tile is in 2 reflected pattern */
        {naddr2 = hashref2();           /* additional moves from 0 pattern */
	   //if ((naddr2 + 1) - addr2 < 0 ) {
	      //printf("addr2 inconsistent\n");
	    //}
        naddr = addr0 + addr1 + naddr2 + addr3;  /* total additional moves */
        if (naddr + g < thresh)                 /* naddr + g + 1 <= thresh */
          if ((nadd == 0) ||                   /* goal state is reached or */
              (search(newblank,blank,g+1,nadd0,add1,add2,add3,addr0,addr1,naddr2,addr3)))
            {path[g] = tile;                          /* record tile moved */
            return (1);}}}                            /* exit with success */

    if (whichpat[tile] == 1)                       /* tile is in 1 pattern */
      {nadd1 = hash1();                 /* additional moves from 1 pattern */
	 //if ((nadd1 + 1) - add1 < 0 ) {
	    // printf("add1 inconsistent\n");
	    //}
      nadd = add0 + nadd1 + add2 + add3;     /* new moves from reg pattern */
      if (nadd + g >= thresh) goto next; /*already > thresh, skip this move*/

      if (whichrefpat[tile] == 2)        /* tile is in 2 reflected pattern */
        {naddr2 = hashref2();           /* additional moves from 0 pattern */
	   //if ((naddr2 + 1) - addr2 < 0 ) {
	      //printf("addr2 inconsistent\n");
	    //}
        naddr = addr0 + addr1 + naddr2 + addr3;  /* total additional moves */
        if (naddr + g < thresh)                 /* naddr + g + 1 <= thresh */
          if ((nadd == 0) ||                   /* goal state is reached or */
              (search(newblank,blank,g+1,add0,nadd1,add2,add3,addr0,addr1,naddr2,addr3)))
            {path[g] = tile;                          /* record tile moved */
            return (1);}}                             /* exit with success */
      else                               /* tile is in 3 reflected pattern */
        {naddr3 = hashref3();           /* additional moves from 0 pattern */
	   //if ((naddr3 + 1) - addr3 < 0 ) {
	      //printf("addr3 inconsistent\n");
	    //}
        naddr = addr0 + addr1 + addr2 + naddr3;  /* total additional moves */
        if (naddr + g < thresh)                 /* naddr + g + 1 <= thresh */
          if ((nadd == 0) ||                   /* goal state is reached or */
              (search(newblank,blank,g+1,add0,nadd1,add2,add3,addr0,addr1,addr2,naddr3)))
            {path[g] = tile;                          /* record tile moved */
            return (1);}}}                            /* exit with success */

    if (whichpat[tile] == 2)                       /* tile is in 2 pattern */
      {nadd2 = hash2();                 /* additional moves from 0 pattern */
	 //if ((nadd2 + 1) - add2 < 0 ) {
	      //printf("add2 inconsistent\n");
	    //}
      nadd = add0 + add1 + nadd2 + add3;     /* new moves from reg pattern */
      if (nadd + g >= thresh) goto next; /*already > thresh, skip this move*/

      if (whichrefpat[tile] == 0)        /* tile is in 0 reflected pattern */
        {naddr0 = hashref0();           /* additional moves from 0 pattern */
	   //if ((naddr0 + 1) - addr0 < 0 ) {
	      //printf("addr0 inconsistent\n");
	    //}
        naddr = naddr0 + addr1 + addr2 + addr3;  /* total additional moves */
        if (naddr + g < thresh)                 /* naddr + g + 1 <= thresh */
          if ((nadd == 0) ||                   /* goal state is reached or */
              (search(newblank,blank,g+1,add0,add1,nadd2,add3,naddr0,addr1,addr2,addr3)))
            {path[g] = tile;                          /* record tile moved */
            return (1);}}                             /* exit with success */
      else                               /* tile is in 1 reflected pattern */
        {naddr1 = hashref1();           /* additional moves from 0 pattern */
	   //if ((naddr1 + 1) - addr1 < 0 ) {
	      //printf("addr1 inconsistent\n");
	    //}
        naddr = addr0 + naddr1 + addr2 + addr3;  /* total additional moves */
        if (naddr + g < thresh)                 /* naddr + g + 1 <= thresh */
          if ((nadd == 0) ||                   /* goal state is reached or */
              (search(newblank,blank,g+1,add0,add1,nadd2,add3,addr0,naddr1,addr2,addr3)))
            {path[g] = tile;                          /* record tile moved */
            return (1);}}}                            /* exit with success */

    if (whichpat[tile] == 3)                       /* tile is in 3 pattern */
      {nadd3 = hash3();                 /* additional moves from 3 pattern */
	 //if ((nadd3 + 1) - add3 < 0 ) {
	      //printf("add3 inconsistent\n");
	    //}
      nadd = add0 + add1 + add2 + nadd3;     /* new moves from reg pattern */
      if (nadd + g >= thresh) goto next; /*already > thresh, skip this move*/

      if (whichrefpat[tile] == 1)        /* tile is in 1 reflected pattern */
        {naddr1 = hashref1();           /* additional moves from 0 pattern */
	   //if ((naddr1 + 1) - addr1 < 0 ) {
	      //printf("addr1 inconsistent\n");
	    //}
        naddr = addr0 + naddr1 + addr2 + addr3;  /* total additional moves */
        if (naddr + g < thresh)                 /* naddr + g + 1 <= thresh */
          if ((nadd == 0) ||                   /* goal state is reached or */
              (search(newblank,blank,g+1,add0,add1,add2,nadd3,addr0,naddr1,addr2,addr3)))
            {path[g] = tile;                          /* record tile moved */
            return (1);}}                             /* exit with success */
      else                               /* tile is in 3 reflected pattern */
        {naddr3 = hashref3();           /* additional moves from 0 pattern */
	   //if ((naddr3 + 1) - addr3 < 0 ) {
	      //printf("addr3 inconsistent\n");
	    //}
        naddr = addr0 + addr1 + addr2 + naddr3;  /* total additional moves */
        if (naddr + g < thresh)                 /* naddr + g + 1 <= thresh */
          if ((nadd == 0) ||                   /* goal state is reached or */
              (search(newblank,blank,g+1,add0,add1,add2,nadd3,addr0,addr1,addr2,naddr3)))
            {path[g] = tile;                          /* record tile moved */
            return (1);}}}                            /* exit with success */

next:s[newblank] = tile;             /* undo current move before doing next */
    inv[tile] = newblank;}
 return (0);}                                          /* exit with failure */

/* READFILE reads in a pattern database file, and expands it into a direct-access database. */
// Modified to take uncompressed database

readfile (table)

unsigned char table[TABLESIZE];      /* direct-access pattern database array */

{
   for (size_t i = 0; i < TABLESIZE; ++i) {
    table[i] = getc(infile);
   }
}

/* MAX takes two integer arguments and returns their maximum. */

max (x, y)

    int x, y;

{if (x >= y) return (x);
 else return (y);}

/* Main program does the initialization, inputs an initial state, solves it,
  and prints the solution. */

main ()

{int success;                          /* boolean flag for success of search */
 int blank;                                     /* initial position of blank */
 int initeval;                        /* manhattan distance of initial state */
 int problem;                                            /* problem instance */
 int index;                                       /* index to tile positions */
 int i;                                                     /* utility index */
 int inith0, inith1, inith2, inith3;      /* initial pattern database values */
 int initrh0, initrh1, initrh2, initrh3; /* initial reflected database values*/
 int hashval;                /* hash value of initial state for each pattern */

 initops ();                                    /* initialize operator table */

 infile = fopen("pat24.1256712.tab", "rb"); /* read 6-tile pattern database */
 readfile (h0);         /* read database and expand into direct-access array */
 fclose(infile);
 printf ("pattern 1 2 5 6 7 12 read in\n");
   

 infile = fopen("pat24.34891314.tab", "rb"); /* read 6-tile pattern database */
 readfile (h1);         /* read database and expand into direct-access array */
 fclose(infile);
 printf ("pattern 3 4 8 9 13 14 read in\n");

 grandtotal = 0;

  blank = input(s);                                 /* input initial state */
                     
  for (index = 0; index < SIZE; index++)          /* for each tile position */
    {inv[s[index]] = index;                     /* initialize inverse state */
    printf ("%2d ", s[index]);}                  /* print out initial state */
  printf ("\n");

  inith0 = hash0();                                  /* new heuristic value */
  inith1 = hash1();                                  /* new heuristic value */
  inith2 = hash2();                                  /* new heuristic value */
  inith3 = hash3();                                  /* new heuristic value */

  initrh0 = hashref0();                              /* new heuristic value */
  initrh1 = hashref1();                              /* new heuristic value */
  initrh2 = hashref2();                              /* new heuristic value */
  initrh3 = hashref3();                              /* new heuristic value */

  thresh = max(inith0+inith1+inith2+inith3,initrh0+initrh1+initrh2+initrh3);
  total = 0;                            /* initialize total nodes generated */

  do                      /* depth-first iterations until solution is found */
    {generated = 0;            /* initialize number generated per iteration */
    success = search(blank, -1, 0, inith0, inith1, inith2, inith3,
                     initrh0, initrh1, initrh2, initrh3);
    printf ("%3d %12.f\n", thresh, generated);           /* nodes per iteration */
    fflush(stdout);
    total = total + generated;     /* keep track of total nodes per problem */
    thresh += 2;}                      /* threshold always increases by two */
  while (!success);                              /* until solution is found */

  grandtotal += total;
  printf ("%d %d %.f\n\n", problem, thresh-2, total);
  for (i = 0; i < thresh-2; i++)                        /* print out solution */
    printf ("%d ", path[i]);
  printf ("\n\n");
  fflush(stdout);

 printf ("grandtotal: %.f\n", grandtotal);}
