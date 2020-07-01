/******************************************************************************
  A template program for developing a GAP solver. Subroutines to read instance
  data and compute the cost of a given solution are included.

  This program can also be used to compute the cost and check the feasibility
  of a solution given from a file. The format of a file is:
  for each job j from 1 to n in this order, the index of the agent (the value
  should be given as values from [1, m]) to which j is assigned. For example,
  if n=4 and m=3, and jobs 1, 2, 3 and 4 are assigned to agents 2, 1, 3 and 1,
  respectively, then the data in the file should be as follows:  2 1 3 1.

  NOTE: Index i of agents ranges from 0 to m-1, and
        index j of jobs   ranges from 0 to n-1 in the program,
	while in the solution file,
	index i of agents ranges from 1 to m, and
        index j of jobs   ranges from 1 to n in the program.
	Sorry for the confusion.

  If you would like to use various parameters, it might be useful to modify
  the definition of struct "Param" and mimic the way the default value of
  "timelim" is given and how its value is input from the command line.
******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "cpu_time.c"

/***** default values of parameters ******************************************/
#define	TIMELIM	300	/* the time limit for the algorithm in seconds */
#define	GIVESOL	0	/* 1: input a solution; 0: do not give a solution */

#define RAND_MAX 0x7fffffff

typedef struct {
  int		timelim;	/* the time limit for the algorithm in secs. */
  int		givesol;	/* give a solution (1) or not (0) */
  /* Never modify the above two lines.  */
  /* You can add more components below. */
} Param;			/* parameters */

typedef struct {
  int	n;	/* number of jobs */
  int	m;	/* number of agents */
  int	**c;	/* cost matrix c_{ij} */
  int	**a;	/* resource requirement matrix a_{ij} */
  int	*b;	/* available amount b_i of resource for each agent i */
} GAPdata;	/* data of the generalized assignment problem */

typedef struct {
  double	timebrid;	/* the time before reading the instance data */
  double	starttime;	/* the time the search started */
  double	endtime;	/* the time the search ended */
  int		*bestsol;	/* the best solution found so far */
  int   *tempsol; 
  int   *tempBestsol;
  int bestcost;
  /* Never modify the above four lines. */
  /* You can add more components below. */
} Vdata;		/* various data often necessary during the search */

/*************************** functions ***************************************/
void copy_parameters(int argc, char *arcv[], Param *param);
void read_instance(GAPdata *gapdata);
void prepare_memory(Vdata *vdata, GAPdata *gapdata);
void free_memory(Vdata *vdata, GAPdata *gapdata);
void read_sol(Vdata *vdata, GAPdata *gapdata);
void recompute_cost(Vdata *vdata, GAPdata *gapdata);
void *malloc_e(size_t size);
int *restB(Vdata *vdata, GAPdata *gapdata);
int isfeasible(Vdata *vdata, GAPdata *gapdata);
int *computeCost(Vdata *vdata, GAPdata *gapdata);

/***** check the feasibility and recompute the cost **************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void recompute_cost(Vdata *vdata, GAPdata *gapdata)
{
  int	i, j;		/* indices of agents and jobs */
  int	*rest_b;	/* the amount of resource available at each agent */
  int	cost, penal;	/* the cost; the penalty = the total capacity excess */
  int	temp;		/* temporary variable */

  rest_b = (int *) malloc_e(gapdata->m * sizeof(int));
  cost = penal = 0;
  for(i=0; i<gapdata->m; i++){rest_b[i] = gapdata->b[i];}
  for(j=0; j<gapdata->n; j++){
    rest_b[vdata->bestsol[j]] -= gapdata->a[vdata->bestsol[j]][j];
    cost += gapdata->c[vdata->bestsol[j]][j];
  }
  for(i=0; i<gapdata->m; i++){
    temp = rest_b[i];
    if(temp<0){penal -= temp;}
  }
  printf("recomputed cost = %d\n", cost);
  if(penal>0){
    printf("INFEASIBLE!!\n");
    printf(" resource left:");
    for(i=0; i<gapdata->m; i++){printf(" %3d", rest_b[i]);}
    printf("\n");
  }
  printf("time for the search:       %7.2f seconds\n",
	 vdata->endtime - vdata->starttime);
  printf("time to read the instance: %7.2f seconds\n",
	 vdata->starttime - vdata->timebrid);

  free((void *) rest_b);
}

int *restB(Vdata *vdata, GAPdata *gapdata){
  int	i, j;		/* indices of agents and jobs */
  int	*rest_b;	/* the amount of resource available at each agent */
  int	cost, penal;	/* the cost; the penalty = the total capacity excess */
  int	temp;		/* temporary variable */
  rest_b = (int *) malloc_e(gapdata->m * sizeof(int));
  cost = penal = 0;
  for(i=0; i<gapdata->m; i++){rest_b[i] = gapdata->b[i];}
  for(j=0; j<gapdata->n; j++){
    rest_b[vdata->bestsol[j]] -= gapdata->a[vdata->bestsol[j]][j];
    cost += gapdata->c[vdata->bestsol[j]][j];
  }
  for(i=0; i<gapdata->m; i++){
    temp = rest_b[i];
    if(temp<0){penal -= temp;}
  }
  return rest_b;
  free((void *) rest_b);
}
int isfeasible(Vdata *vdata, GAPdata *gapdata){
  // 1ならok 0ならだめ
  int	i, j;		/* indices of agents and jobs */
  int	*rest_b;	/* the amount of resource available at each agent */
  int	cost, penal;	/* the cost; the penalty = the total capacity excess */
  int	temp;		/* temporary variable */
  rest_b = (int *) malloc_e(gapdata->m * sizeof(int));
  cost = penal = 0;
  for(i=0; i<gapdata->m; i++){rest_b[i] = gapdata->b[i];}
  for(j=0; j<gapdata->n; j++){
    rest_b[vdata->bestsol[j]] -= gapdata->a[vdata->bestsol[j]][j];
    cost += gapdata->c[vdata->bestsol[j]][j];
  }
  for(i=0; i<gapdata->m; i++){
    temp = rest_b[i];
    if(temp<0){penal -= temp;}
  }
  if(penal>0){
    return 0;
  }else{
    return 1;
  }
  
  free((void *) rest_b);
}

int *computeCost(Vdata *vdata, GAPdata *gapdata){
  // 1番目にfeasibleかどうか，2番目にcostを返す
  int	i, j;		/* indices of agents and jobs */
  int	*rest_b;	/* the amount of resource available at each agent */
  int	cost, penal;	/* the cost; the penalty = the total capacity excess */
  int	temp;		/* temporary variable */
  int *ReturnArray;
  rest_b = (int *) malloc_e(gapdata->m * sizeof(int));
  cost = penal = 0;

  ReturnArray = (int *) malloc_e(2);
  for(i=0; i<gapdata->m; i++){rest_b[i] = gapdata->b[i];}
  for(j=0; j<gapdata->n; j++){
    rest_b[vdata->tempsol[j]] -= gapdata->a[vdata->tempsol[j]][j];
    cost += gapdata->c[vdata->tempsol[j]][j];
  }
  for(i=0; i<gapdata->m; i++){
    temp = rest_b[i];
    if(temp<0){penal -= temp;}
  }
  ReturnArray[0] = penal;
  ReturnArray[1] = cost;
  return ReturnArray;
  free((void *) rest_b);
  free((void *) ReturnArray);

}

/***** read a solution from STDIN ********************************************/
void read_sol(Vdata *vdata, GAPdata *gapdata)
{
  int	j;		/* index of jobs */
  int	value_read;	/* the value read by fscanf */
  FILE	*fp=stdin;	/* set fp to the standard input */

  for(j=0; j<gapdata->n; j++){
    fscanf(fp, "%d", &value_read);
    /* change the range of agents from [1, m] to [0, m-1] */
    vdata->bestsol[j] = value_read - 1;
  }
}

/***** prepare memory space **************************************************/
/***** Feel free to modify this subroutine. **********************************/
void prepare_memory(Vdata *vdata, GAPdata *gapdata)
{
  int j;

  vdata->bestsol = (int *)  malloc_e(gapdata->n * sizeof(int));
  vdata->tempsol = (int *)  malloc_e(gapdata->n * sizeof(int));
  vdata->tempBestsol = (int *)  malloc_e(gapdata->n * sizeof(int));
  /* the next line is just to avoid confusion */
  for(j=0; j<gapdata->n; j++){vdata->bestsol[j] = 0;}
  for(j=0; j<gapdata->n; j++){vdata->tempsol[j] = 0;}
  for(j=0; j<gapdata->n; j++){vdata->tempBestsol[j] = 0;}
}

/***** free memory space *****************************************************/
/***** Feel free to modify this subroutine. **********************************/
void free_memory(Vdata *vdata, GAPdata *gapdata)
{
  free((void *) vdata->bestsol);
  free((void *) vdata->tempsol);
  free((void *) vdata->tempBestsol);
  free((void *) gapdata->c[0]);
  free((void *) gapdata->c);
  free((void *) gapdata->a[0]);
  free((void *) gapdata->a);
  free((void *) gapdata->b);
}

/***** read the instance data ************************************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void read_instance(GAPdata *gapdata)
{
  int	i, j;		/* indices of agents and jobs */
  int	value_read;	/* the value read by fscanf */
  FILE	*fp=stdin;	/* set fp to the standard input */

  /* read the number of agents and jobs */
  fscanf(fp, "%d", &value_read);	/* number of agents */
  gapdata->m = value_read;
  fscanf(fp,"%d",&value_read);		/* number of jobs */
  gapdata->n = value_read;

  /* initialize memory */
  gapdata->c    = (int **) malloc_e(gapdata->m * sizeof(int *));
  gapdata->c[0] = (int *)  malloc_e(gapdata->m * gapdata->n * sizeof(int));
  for(i=1; i<gapdata->m; i++){gapdata->c[i] = gapdata->c[i-1] + gapdata->n;}
  gapdata->a    = (int **) malloc_e(gapdata->m * sizeof(int *));
  gapdata->a[0] = (int *)  malloc_e(gapdata->m * gapdata->n * sizeof(int));
  for(i=1; i<gapdata->m; i++){gapdata->a[i] = gapdata->a[i-1] + gapdata->n;}
  gapdata->b    = (int *)  malloc_e(gapdata->m * sizeof(int));

  /* read the cost coefficients */   
  for(i=0; i<gapdata->m; i++){    
    for(j=0; j<gapdata->n; j++){
      fscanf(fp, "%d", &value_read);
      gapdata->c[i][j] = value_read;
    }
  }

  /* read the resource consumption */
  for(i=0; i<gapdata->m; i++){
    for(j=0; j<gapdata->n; j++){
      fscanf(fp, "%d", &value_read);
      gapdata->a[i][j] = value_read;
    }
  }

  /* read the resource capacity */
  for(i=0; i<gapdata->m; i++){    
    fscanf(fp,"%d", &value_read);
    gapdata->b[i] = value_read;
  }
}

/***** copy and read the parameters ******************************************/
/***** Feel free to modify this subroutine. **********************************/
void copy_parameters(int argc, char *argv[], Param *param)
{
  int i;

  /**** copy the parameters ****/
  param->timelim = TIMELIM;
  param->givesol = GIVESOL;
  /**** read the parameters ****/
  if(argc>0 && (argc % 2)==0){
    printf("USAGE: ./gap [param_name, param_value] [name, value]...\n");
    exit(EXIT_FAILURE);}
  else{
    for(i=1; i<argc; i+=2){
      if(strcmp(argv[i],"timelim")==0) param->timelim = atoi(argv[i+1]);
      if(strcmp(argv[i],"givesol")==0) param->givesol = atoi(argv[i+1]);
    }
  }
}

/***** malloc with error check ***********************************************/
void *malloc_e( size_t size ) {
  void *s;
  if ( (s=malloc(size)) == NULL ) {
    fprintf( stderr, "malloc : Not enough memory.\n" );
    exit( EXIT_FAILURE );
  }
  return s;
}



void ShiftShuffle(int arr[],int n){
  for (int i=n;i>1;i--){
    int a=i-1;
    int b=rand()%i;
    int tmp = arr[b];
    arr[b] = arr[a];
    arr[a] = tmp;
  }
  return;
}
void SwapShuffle(int arr[],int n){
  for (int i=n;i>1;i--){
    int a=i-1;
    int b=rand()%i;
    int tmp = arr[b];
    arr[b] = arr[a];
    arr[a] = tmp;
  }
  return;
}

int accept_prob(double delta, double T){
  if (delta<0){
    return 1;
  }else if ((double)exp(-delta/T) > (double)rand()/RAND_MAX){
    return 1;
  }else{
    return 0;
  }
}


/***** main ******************************************************************/
int main(int argc, char *argv[])
{
  Param		param;		/* parameters */
  GAPdata	gapdata;	/* GAP instance data */
  Vdata		vdata;		/* various data often needed during search */

  vdata.timebrid = cpu_time();
  copy_parameters(argc, argv, &param);
  read_instance(&gapdata);
  prepare_memory(&vdata, &gapdata);
  if(param.givesol==1){read_sol(&vdata, &gapdata);}
  vdata.starttime = cpu_time();

  /*
    Write your program here. Of course you can add your subroutines
    outside main(). At this point, the instance data is stored in "gapdata".
	gapdata.n	number of jobs n
	gapdata.m	number of agents m
	gapdata.c[i][j]	cost c_{ij} 
	gapdata.a[i][j]	resource requirement a_{ij} 
	gapdata.b[i]	available amount b_i of resource at agent i
    Note that i ranges from 0 to m-1, and j ranges from 0 to n-1. Note also
    that  you should write, e.g., "gapdata->c[i][j]" in your subroutines.
    Store your best solution in vdata.bestsol, then "recompute_cost" will
    compute its cost and its feasibility. The format of vdata.bestsol is:
    For each job j from 0 to n-1 in this order, the index of the agent 
    (the value should be given as values from [0, m-1]) to which j is
    assigned. For example, if n=4 and m=3, and jobs 0, 1, 2 and 3 are
    assigned to agents 1, 0, 2 and 0, respectively, then vdata.bestsol
    should be as follows:  
	vdata.bestsol[0] = 1
	vdata.bestsol[1] = 0
	vdata.bestsol[2] = 2
	vdata.bestsol[3] = 0.
    Note that you should write "vdata->bestsol[j]" in your subroutines.
  */
  
  int b[gapdata.m];
  int a[gapdata.m][gapdata.n];
  for(int i=0;i<gapdata.n;i++){
    for (int j=0;j<gapdata.m;j++){
      a[j][i] = gapdata.a[j][i];
    }
  }
  for (int i=0;i<gapdata.m;i++){
    b[i] = gapdata.b[i];
  }
  int F[gapdata.n][gapdata.m];
  for (int iternum=0;iternum<gapdata.n;iternum++){
    for (int i=0;i<gapdata.n;i++){
      for (int j=0;j<gapdata.m;j++){
        if (a[j][i] <= b[j]){
          F[i][j] = j;
        }else{
          F[i][j] = -1;
        }
      }
    }
    int MaxDiff = 0;
    int MaxDiffIdx = 0;
    int DeleteIdx = 0;
    for (int i=0;i<gapdata.n;i++){
      int MinIdx = 0;
      int MinVal = 99999999;
      int SecondMinVal=99999999;
      int MaxVal = 0;
      for (int j=0;j<gapdata.m;j++){
        if (F[i][j] != -1){
          if (a[j][i] < MinVal){
            MinVal = a[j][i];
            MinIdx = j;
          }
          else if (a[j][i] < SecondMinVal){
            MaxVal = a[j][i];
            SecondMinVal = a[j][i];
          } 
        }
      }
      if (MinVal != 99999999 && SecondMinVal == 99999999){
        MaxDiff = 0;
        MaxDiffIdx = i;
        DeleteIdx = MinIdx;
      }else if (MinVal != 99999999 && SecondMinVal != 99999999){
        if (SecondMinVal-MinVal > MaxDiff){
          MaxDiff = SecondMinVal-MinVal;
          MaxDiffIdx = i;
          DeleteIdx = MinIdx;
        }
      }
    }
    if (a[DeleteIdx][MaxDiffIdx] == 999){
      int Maxb = b[0];
      for (int i=1;i<gapdata.m;i++){
        if (Maxb < b[i]) {
          DeleteIdx = i;
          Maxb = b[i];
        }
      }

    }
    vdata.bestsol[MaxDiffIdx] = DeleteIdx;
    b[DeleteIdx] -= a[DeleteIdx][MaxDiffIdx];
    for (int i=0;i<gapdata.m;i++){
      a[i][MaxDiffIdx] = 999;
    }

    


  }
  int ok = isfeasible(&vdata, &gapdata);
  while (ok == 0){
    int MaxIdx=-1;
    int MaxRest = -999999;
    int MinIdx=-1;
    int MinRest = 999999;
    int* restArray = restB(&vdata, &gapdata);
    for (int i=0;i<gapdata.m;i++){
      if (restArray[i] < MinRest){
        MinRest=restArray[i];
        MinIdx = i;
      }
      if (restArray[i] > MaxRest){
        MaxRest = restArray[i];
        MaxIdx = i;
      }
    }  
    for (int i=0;i<gapdata.n;i++){
      if (vdata.bestsol[i] == MinIdx) {
        vdata.bestsol[i] = MaxIdx;
        break;
      }
    }
    ok = isfeasible(&vdata, &gapdata);
  }

  // 初期解生成 終了
  for (int i=0;i<gapdata.n;i++){
    vdata.tempsol[i] = vdata.bestsol[i];
  }
  int BestCost = computeCost(&vdata, &gapdata)[1];





  // 局所探索するリストを作成する準備
  int ShiftNeighbor[gapdata.n*gapdata.m][2];
  int SwapNeighbor[gapdata.n*gapdata.n][2];
  int testIdx=0;
  for (int i=0;i<gapdata.m;i++){
    for (int j=0;j<gapdata.n;j++){
      ShiftNeighbor[testIdx][0] = j;
      ShiftNeighbor[testIdx][1] = i;
      testIdx++;
    }
  }
  testIdx=0;
  for (int i=0;i<gapdata.n;i++){
    for (int j=0;j<gapdata.n;j++){
      SwapNeighbor[testIdx][0] = i;
      SwapNeighbor[testIdx][1] = j; 
      testIdx++;
    }
  }
  int ShiftRandomIdx[gapdata.n*gapdata.m];
  for (int i=0;i<gapdata.n*gapdata.m;i++){
    ShiftRandomIdx[i] = i;
  }
  int SwapRandomIdx[gapdata.n*gapdata.n];
  for (int i=0;i<gapdata.n*gapdata.n;i++){
    SwapRandomIdx[i] = i;
  }



  int AgentIdx,JobIdx;
  int FirstAgentIdx,SecondAgentIdx,FirstJobIdx,SecondJobIdx;

  int *Result;


  
  /*
  // 局所探索
  int noImprove = 0;
  while(noImprove == 0){
    noImprove = 1;
    ShiftShuffle(ShiftRandomIdx,gapdata.n*gapdata.m);
    for (int i=0;i<gapdata.n*gapdata.m;i++){
      AgentIdx = ShiftNeighbor[ShiftRandomIdx[i]][0];
      JobIdx = ShiftNeighbor[ShiftRandomIdx[i]][1];
      // printf("%d,%d\n",AgentIdx,JobIdx);
      vdata.tempsol[AgentIdx] = JobIdx;
      Result = (int *)malloc(2 * sizeof(int));
      Result = computeCost(&vdata, &gapdata);
      // infeasibleの解を訪れない場合
      if (Result[0] == 0 && Result[1]<BestCost){
        vdata.bestsol[AgentIdx] = vdata.tempsol[AgentIdx];
        BestCost = Result[1];
        noImprove=0;
        printf("%d\n",BestCost);
      }else{
        vdata.tempsol[AgentIdx] = vdata.bestsol[AgentIdx];
      }
      free(Result);  
    }
    printf("shift近傍探索\n");
    SwapShuffle(SwapRandomIdx,gapdata.n*gapdata.n);
    for (int i=0;i<gapdata.n*gapdata.n;i++){
      FirstAgentIdx = SwapNeighbor[SwapRandomIdx[i]][0];
      FirstJobIdx = vdata.tempsol[FirstAgentIdx];
      SecondAgentIdx = SwapNeighbor[SwapRandomIdx[i]][1];
      SecondJobIdx = vdata.tempsol[SecondAgentIdx];
      // printf("%d,%d\n",AgentIdx,JobIdx);
      vdata.tempsol[SecondAgentIdx] = FirstJobIdx;
      vdata.tempsol[FirstAgentIdx] = SecondJobIdx;

      Result = (int *)malloc(2 * sizeof(int));
      Result = computeCost(&vdata, &gapdata);
      if (Result[0] == 0 && Result[1]<BestCost){
        vdata.bestsol[FirstAgentIdx] = vdata.tempsol[FirstAgentIdx];
        vdata.bestsol[SecondAgentIdx] = vdata.tempsol[SecondAgentIdx];
        BestCost = Result[1];
        noImprove=0;
        printf("%d\n",BestCost);
      }else{
        vdata.tempsol[FirstAgentIdx] = vdata.bestsol[FirstAgentIdx];
        vdata.tempsol[SecondAgentIdx] = vdata.bestsol[SecondAgentIdx];
      }
      free(Result);  
    }
    printf("swap近傍探索\n");
    if (noImprove == 1){
      printf("改善なし\n");
    }
  }
  */


  // アニーリング法の初期温度設定
  double T=10;
  int tmp_delta,accepted_cnt,all_cnt;
  int tmp;
  do{
    accepted_cnt = 0;
    all_cnt=0;
    T = T*2;
    while(all_cnt<100){
      AgentIdx = rand() % gapdata.n;
      JobIdx = rand() % gapdata.m;
      while(JobIdx == vdata.tempsol[AgentIdx]){
        JobIdx = rand() % gapdata.m;
      }
      vdata.tempsol[AgentIdx] = JobIdx;
      Result = (int *)malloc(2 * sizeof(int));
      Result = computeCost(&vdata, &gapdata);
      if (Result[0]== 0){
        all_cnt += 1;
        tmp_delta = Result[1]-BestCost;
        if (accept_prob(tmp_delta,T)==1){
          accepted_cnt += 1;
        }
      }
      free(Result);
      vdata.tempsol[AgentIdx] = vdata.bestsol[AgentIdx];
    }
    while(all_cnt<200){
      FirstAgentIdx = rand() % gapdata.n;
      SecondAgentIdx = rand() % gapdata.n;
      while(FirstAgentIdx == SecondAgentIdx){
        SecondAgentIdx = rand() % gapdata.n;
      }
      tmp = vdata.tempsol[FirstAgentIdx];
      vdata.tempsol[FirstAgentIdx] = vdata.tempsol[SecondAgentIdx];
      vdata.tempsol[SecondAgentIdx] = tmp;
      Result = (int *)malloc(2 * sizeof(int));
      Result = computeCost(&vdata, &gapdata);
      if (Result[0]== 0){
        all_cnt += 1;
        tmp_delta = Result[1]-BestCost;
        if (accept_prob(tmp_delta,T)==1){
          accepted_cnt += 1;
        }
      }
      free(Result);
      vdata.tempsol[FirstAgentIdx] = vdata.bestsol[FirstAgentIdx];
      vdata.tempsol[SecondAgentIdx] = vdata.bestsol[SecondAgentIdx];
    }
  }while(accepted_cnt <= 180);

  
  int tempBestCost = BestCost;
  for (int i=0;i<gapdata.n;i++){
    vdata.tempBestsol[i] = vdata.bestsol[i];
  }

  int delta;
  int consecutive_count;
  while(T>1.0){
    consecutive_count = 0;
    // shift局所探索
    ShiftShuffle(ShiftRandomIdx,gapdata.n*gapdata.m);
    for (int i=0;i<gapdata.n*gapdata.m;i++){
      if (consecutive_count >=50){
        break;
      }
      AgentIdx = ShiftNeighbor[ShiftRandomIdx[i]][0];
      JobIdx = ShiftNeighbor[ShiftRandomIdx[i]][1];
      vdata.tempsol[AgentIdx] = JobIdx;
      Result = (int *)malloc(2 * sizeof(int));
      Result = computeCost(&vdata, &gapdata);
      if (Result[0] == 0){
        delta = Result[1]-tempBestCost;
        if (accept_prob(delta,T)==1){
          consecutive_count += 1;
          vdata.tempBestsol[AgentIdx] = vdata.tempsol[AgentIdx];
          tempBestCost = Result[1];
          if (delta<0){
            for (int i=0;i<gapdata.n;i++){
              vdata.bestsol[i] = vdata.tempBestsol[i];
            }
          }
        }else{
          vdata.tempsol[AgentIdx] = vdata.tempBestsol[AgentIdx];
          consecutive_count = 0;
        }
      }else{
        vdata.tempsol[AgentIdx] = vdata.tempBestsol[AgentIdx];
      }
      free(Result);  
    }

    consecutive_count = 0;

    // swap近傍
    SwapShuffle(SwapRandomIdx,gapdata.n*gapdata.n);
    for (int i=0;i<gapdata.n*gapdata.n;i++){
      if (consecutive_count >=50){
        break;
      }
      FirstAgentIdx = SwapNeighbor[SwapRandomIdx[i]][0];
      FirstJobIdx = vdata.tempsol[FirstAgentIdx];
      SecondAgentIdx = SwapNeighbor[SwapRandomIdx[i]][1];
      SecondJobIdx = vdata.tempsol[SecondAgentIdx];
      vdata.tempsol[SecondAgentIdx] = FirstJobIdx;
      vdata.tempsol[FirstAgentIdx] = SecondJobIdx;
      Result = (int *)malloc(2 * sizeof(int));
      Result = computeCost(&vdata, &gapdata);
      if (Result[0] == 0){
        delta = Result[1]-tempBestCost;
        if (accept_prob(delta,T)==1){
          consecutive_count += 1;
          vdata.tempBestsol[FirstAgentIdx] = vdata.tempsol[FirstAgentIdx];
          vdata.tempBestsol[SecondAgentIdx] = vdata.tempsol[SecondAgentIdx];
          tempBestCost = Result[1];
          if (delta<0){
            for (int i=0;i<gapdata.n;i++){
              vdata.bestsol[i] = vdata.tempsol[i];
            }
          }
        }else{
          consecutive_count = 0;
          vdata.tempsol[FirstAgentIdx] = vdata.tempBestsol[FirstAgentIdx];
          vdata.tempsol[SecondAgentIdx] = vdata.tempBestsol[SecondAgentIdx];
        }
      }else{
        vdata.tempsol[FirstAgentIdx] = vdata.tempBestsol[FirstAgentIdx];
        vdata.tempsol[SecondAgentIdx] = vdata.tempBestsol[SecondAgentIdx];
      }
      free(Result);  
    }
    if ((double)cpu_time() - vdata.starttime  > TIMELIM){
      break;
    }


    T= T*0.97;
    // printf("T %f\n ",T);
  }




  // for (int i=0;i<gapdata.n;i++){
  //   printf("%d ",vdata.bestsol[i]);
  // }
  // printf("\n");
  vdata.endtime = cpu_time();
  recompute_cost(&vdata, &gapdata);
  free_memory(&vdata, &gapdata);

  return EXIT_SUCCESS;
}
