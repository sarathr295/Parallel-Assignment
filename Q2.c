 #include <iostream>
   #include <algorithm>
   #include <math.h>
   #include <iomanip>
   #include <mpi.h>
  using namespace std;
  
   const int N=8000;
   double Ab[N][N+1];
  double oriAb[N][N+1]; //original Matrix Ab for testing the solution
  double x[N]={0}; //the solution vector
 
  void PrintMatrix();
  void SetMatrix();
  void TestSolution();
  
  int main(int argc, char **argv){
    int id, nproc, id_from, id_to;
    MPI_Status status;
 
  MPI_Init(&argc, &argv);
     MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &nproc);
   
    int irow,jrow, i,j, jmax;
    double t, amul, item, item_max;
    time_t ts, te;
  
    if(id==0){
    SetMatrix();
    cout<<"The original matrix:"<<endl;
   PrintMatrix();
  
    ts=MPI_Wtime();
    }
   MPI_Bcast(Ab, N*(N+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
   for(irow=0;irow<N-1;++irow){
    id_from=(N-1-irow)%nproc;
 
      jmax=irow;
      item_max=fabs(Ab[irow][irow]);
      MPI_Bcast(&item_max, 1, MPI_DOUBLE, id_from, MPI_COMM_WORLD);
 
     for(jrow=N-1-id;jrow>=irow+1;jrow-=nproc)
       if(fabs(Ab[jrow][irow])>item_max){
  	jmax=jrow;
 	item_max=fabs(Ab[jrow][irow]);
      }
      
     struct{
       double item;
       int row;
     } p, tmp_p;
     p.item=item_max; p.row=jmax;
     MPI_Allreduce(&p, &tmp_p, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
     jmax=tmp_p.row;
     id_to=(N-1-jmax)%nproc;
 
     if(id_from==id_to){ //replacement at the same slave: no communication
       if(id==id_from && jmax!=irow)
 	for(j=0;j<N+1;++j){
 	  item=Ab[irow][j];
 	  Ab[irow][j]=Ab[jmax][j];
 	  Ab[jmax][j]=item;
 	}
     }else{ //communication is needed for pivot exchange
       if(id==id_to){
 	for(j=irow;j<N+1;++j)
 	  Ab[irow][j]=Ab[jmax][j];
 	MPI_Send(&Ab[irow], N+1, MPI_DOUBLE, id_from, 77, MPI_COMM_WORLD);
 	MPI_Recv(&Ab[jmax], N+1, MPI_DOUBLE, id_from, 88, 
 		 MPI_COMM_WORLD, &status);
 	}
       if(id==id_from){
  	for(j=irow;j<N+1;++j)
 	  Ab[jmax][j]=Ab[irow][j];
 	MPI_Recv(&Ab[irow], N+1, MPI_DOUBLE, id_to, 77, 
 		 MPI_COMM_WORLD, &status);
 	MPI_Send(&Ab[jmax], N+1, MPI_DOUBLE, id_to, 88, MPI_COMM_WORLD);
       }
     }
 
    MPI_Bcast(&Ab[irow], N+1, MPI_DOUBLE, id_from, MPI_COMM_WORLD);
    for(jrow=N-1-id;jrow>=irow+1;jrow-=nproc){
       amul=Ab[jrow][irow] * t;

      for(j=irow;j<N+1;++j)
	Ab[jrow][j] += amul * Ab[irow][j];
  }
  }

 if(id==0){
   cout<<"The upper triangular matrix:"<<endl;
   PrintMatrix(); 
    for(irow=N-1;irow>=0;--irow){
     x[irow]= - Ab[irow][N]/Ab[irow][irow];
     for(jrow=0;jrow<irow;++jrow){
	Ab[jrow][N] += x[irow] * Ab[jrow][irow];
	Ab[jrow][irow]=0;
      }
    }
    te=MPI_Wtime();
    cout<<"time elapsed: "<<te-ts<<endl;
    
     cout<<"The solution matrix:"<<endl;
    PrintMatrix();
    TestSolution();
  }


   MPI_Finalize();
 }
 void SetMatrix(){
  int i,j;
   srand(time(NULL));
   for(i=0;i<N;++i)
    for(j=0;j<N+1;++j)
      oriAb[i][j]=Ab[i][j]=rand()%10;
 }
 
 void PrintMatrix(){

  int i,j;
   if(N>20){cout<<"Too big to display!"<<endl;return;}
  cout.precision(1);
   for(i=0;i<N+1;++i)cout<<"------";
  cout<<fixed<<"-----"<<endl;
   for(i=0;i<N;++i){
    cout<<"| ";
    for(j=0;j<N;++j)
      cout<<setw(5)<<Ab[i][j]<<" ";
    cout<<"| "<<setw(5)<<Ab[i][j];
    cout<<" |  x["<<setw(2)<<i<<"] = "<<setw(5)<<x[i]<<endl;
   }
   for(i=0;i<N+1;++i)cout<<"------";
   cout<<"-----"<<endl;
 }
 
 void TestSolution(){
   int i,j;
   double diff, sum;
   cout.precision(20);
  for(i=0;i<N;++i){
    sum=0;
    for(j=0;j<N;++j)
      sum += x[j] * oriAb[i][j];
    diff=sum+oriAb[i][N];
    if(diff>0.0001 || diff<-0.0001)
     cout<<"ERROR! "<<sum<<" ~ "<<oriAb[i][N]<<", diff:"<<diff<<endl;
     if(N<50){
      cout<<setw(4)<<sum<<" ~ "<<setw(4)<<oriAb[i][N];
       cout<<", diff:"<<setw(4)<<fixed<<diff<<endl;
    }
 }
 }
