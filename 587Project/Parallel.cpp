/*
Zaina Connect 4 Application to be run in Parallel across Processors Using MPI & C++

At each given state, the Processors should look through all the game trees and decide what would be the best move for player 1.
Alpha Beta pruning & min-max algorithm to ensure the entire tree doesn't need to be traversed while searching for an optimal path.
Load Balancing needs to be done in order to ensure all the processors have some amount of work to be done at every point in time.
*/

/*
Functions:
1. stateGenerator() xx
3. checkDraw() xx
4. checkForWin() xx
5. displayBoard() xx
6. evalBoard() xx 
7. generateGlobalQueue() xx 
	this can choose to have atleast double OR triple the number of 	startStates for the workers to take from the manager 


7. checkIfVisited() //takes in the state's matrix and returns true if this has been visited. // if it was visited, also makes sure the path's are also copied, is this even required?

The manager has a list of 3 level down states
	check the eval value before sending it to the recursive function & before adding it to the queue as well
	if the eval in any level where the player 'X' just played is 999
		If its the first level itself that is the move to be given out immediately
		If its at a deeper level, we need to compare against the B value of the previous level, and return that state for which the B value was the lowest
*/

#include <mpi.h>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <tuple>
#include <vector>
#include <utility> 
#include <functional>
#include <string>
#include <unordered_map>
#include <queue>
#include <deque>
using namespace std;

const int numrows = 6;
const int numcols = 7;
const int evaluationTable[6][7] = {{3, 4, 5, 7, 5, 4, 3}, {4, 6, 8, 10, 8, 6, 4}, {5, 8, 11, 13, 11, 8, 5}, {5, 8, 11, 13, 11, 8, 5}, {4, 6, 8, 10, 8, 6, 4}, {3, 4, 5, 7, 5, 4, 3}};

double starttime=0;
double endtime=0;
double timeTaken=0;

std::string level;
int pid, numProcs; //processor number & total no of processors
int globalBestVal;

std::unordered_map<std::string,bool> globalHashMap; //Maps to the String of the matrix state & the Bool that gives its visited value
std::unordered_map<std::string,bool> localHashMap;

struct state {

    int path [12];
    int board [numrows][numcols];
    int player;
    int boardSum;
    int colHeight [numcols];
	//write a giveHash() function inside this, that takes the matrix, and sets a hash value for it in the global hash map,
};

#define NBLOCKS 5

struct state startState;
struct state localBestState;
struct state globalBestState;

std::deque<struct state> globalQueue;
std::deque<struct state> localQueue;

//gives the hash value for that state
std::string giveHash(struct state temp)
{
	//std::vector<std::vector<int>> board = temp.board;
	std::string hashVal;
    for(int c=0;c<numcols;c++)
    {
		for(int r=0;r<numrows;r++)
		{
			if (temp.board[r][c]==0) //End of elements in that column
				break;
            else if (temp.board[r][c]==1)
            	hashVal+='1';
            else
            	hashVal+='2';
		}
		hashVal+="|";
	}		
    return hashVal;
}


//displays the board 'temp'
void displayBoard(struct state temp)
{
	//display board heading
	cout<<"\n\n|";
	for(int col=0; col<numcols; col++)
		cout<<" "<<col+1<<"  ";
	cout<<"|\n";

	//std::vector<std::vector<int>> board = temp.board;

	for(int row = numrows-1; row>=0; row--)
	{	cout<<"|";
		for(int col=0; col<numcols; col++)
		{
			if(temp.board[row][col] == 1)
				cout<<" X  ";
			else if (temp.board[row][col] == 2)
				cout<<" O  ";
			else
				cout<<" .  ";
		}
		cout<<"|\n";
	}
	cout<<"'";
	for(int col=0; col<numcols; col++)
		cout<<"----";
	cout<<"'\n\n";

	/*Printing the column heights
	cout<<"Column Height : [";
	for(int i =0; i<temp.colHeight.size(); i++)
		cout<<temp.colHeight[i]<<" ";
	cout<<"]"<<endl;*/

	//Printing the path vector
	/*if(temp.path.size())
	{
		cout<<"\n Path Taken : ";
		for (int i=0;i<temp.path.size(); i++)
		{
			//each path[i] is the first path taken which is a pair of pairs
			cout<<"\n\t Player "<<temp.path[i].second<<" : at ["<<temp.path[i].first.first<<","<<temp.path[i].first.second<<"].";
		}
		cout<<endl;
	}*/
}


int initialiseState()
{
	startState.player = 1;
	int sum;
	//allocate memory for our state
	//std::vector<std::vector<int>> board;

	for(int r=0;r<numrows;r++)
		for(int c=0;c<numcols;c++)
			startState.board[r][c]=0;

	if (level.compare("full") == 0)
	{
		startState.board[0][0] = startState.board[0][2] = startState.board[0][3] = startState.board[0][6] = startState.board[1][1] = startState.board[1][2] 
		= startState.board[1][4] = startState.board[1][5] = startState.board[2][0] = startState.board[2][4] = startState.board[2][6] = startState.board[3][1] 
		= startState.board[3][5] = startState.board[4][3] = startState.board[4][6] = 1;
		
		startState.board[0][1] = startState.board[0][4] = startState.board[0][5] = startState.board[1][0] = startState.board[1][3] = startState.board[1][6] 
		= startState.board[2][1] = startState.board[2][2] = startState.board[2][3] = startState.board[2][5] = startState.board[3][2] = startState.board[3][3] 
		= startState.board[3][4] = startState.board[3][6] = startState.board[4][5] = 2;

		startState.boardSum = 45;
	}

	else if (level.compare("middle") == 0)
	{
		startState.board[0][0] = startState.board[0][1] = startState.board[0][2] = startState.board[0][4] = startState.board[1][2] = startState.board[2][3] = startState.board[2][4] = 1;
		startState.board[0][3] = startState.board[1][0] = startState.board[1][1] = startState.board[1][3] = startState.board[1][4] = startState.board[2][2] = startState.board[3][4] = 2;
		startState.boardSum = 21;
	}
	else if (level.compare("begin") == 0)
	{
		startState.board[0][1] = startState.board[0][3] = startState.board[0][4] = 1;
		startState.board[0][2] = startState.board[0][5] = startState.board[0][6] = 2;
		startState.boardSum = 9;
	}
	else //incorrect Input
	{
		cout<<"\n Incorrect Input, Please start again!";
		return -1;
	}

	
	//set startState.colHeight
	for(int c=0;c<numcols;c++)
	{
		int height = 0;
		for(int r=0;r<numrows;r++)
			if(startState.board[r][c] != 0)
				++height;

		startState.colHeight[c] = height;
		//cout<<"\n Height of column '"<<c<<"' is '"<<height<<endl;
	}

	//startState.board = board;
	//displayBoard(startState);
	return 0;
}


int checkWin(struct state temp)
{
	//std::vector<std::vector<int>> board = temp.board;
	int winner = 0;
	for(int row = 0; row < numrows; row ++)
	{
		for(int col = 0; col < numcols; col ++)
		{
			if(temp.board[row][col] != 0)
			{
				// Across
				if(col < numcols-3)
				{
					if (temp.board[row][col] == temp.board[row][col+1] &&
						temp.board[row][col] == temp.board[row][col+2] &&
						temp.board[row][col] == temp.board[row][col+3])
						{
							winner = temp.board[row][col];
							break;
						}
				}

				// Upwards/downwards
				if(row < numrows-3)
				{
					if (temp.board[row][col] == temp.board[row+1][col] &&
						temp.board[row][col] == temp.board[row+2][col] &&
						temp.board[row][col] == temp.board[row+3][col])
					{
						winner = temp.board[row][col];
						break;
					}
				}

				//Diagonal down-right or up-left
				if(row < numrows-3 && col < numcols-3)
				{
					if (temp.board[row][col] == temp.board[row+1][col+1] &&
						temp.board[row][col] == temp.board[row+2][col+2] &&
						temp.board[row][col] == temp.board[row+3][col+3])
					{
						winner = temp.board[row][col];
						break;
					}
				}
				
				//Diagonal down-left or up-right
				if(row > 3 && col < numcols-3)
				{
					if (temp.board[row][col] == temp.board[row-1][col+1] &&
						temp.board[row][col] == temp.board[row-2][col+2] &&
						temp.board[row][col] == temp.board[row-3][col+3])
					{
						winner = temp.board[row][col];
						break;
					}
				}
			}
		}
		if(winner)
			break;
	}
	return (winner);
}

// int testWin ()
// {
// 	std::vector<std::vector<int>> board = startState.board;
// 	/*board[4][3] = 2;
// 	board[4][5] = 1;
// 	board[1][0] = 1;
// 	board[2][0] = 2;
// 	board[2][4] = 2;
// 	board[3][4] = 1;*/
// 	int winner = checkWin(startState);
// 	if(winner == 0)
// 		cout<<"No Winner!"<<endl;
// 	else
// 		cout<<"\nWinner = Player "<<winner<<". (1 = X & 2 = O)"<<endl;
	
// 	displayBoard(startState);
// 	return 1;
// }

int checkDraw(struct state temp)
{
	if((!checkWin(temp)) && temp.boardSum == 63) //when the temp.boardSum == 63 it means the board is filled with 21 1's and 21 2's = 63.
		return 1;
	return 0;
}

int evalBoard(struct state temp)
{
	//std::vector<std::vector<int>> board = temp.board;
	int winner = checkWin(temp);
	if (winner == 1)
		return +999;
	else if (winner == 2)
		return -999;

	else //No player has won
	{
		if (checkDraw(temp))
			return 0;
		else// AND no draw
		{
			int evalSum = 0;
			for(int r=0;r<numrows;r++)
				for(int c=0;c<numcols;c++)
					if(temp.board[r][c] == 1)
						evalSum += evaluationTable[r][c];
					else if (temp.board[r][c]== 2)
						evalSum -= evaluationTable[r][c];
			return evalSum;
		}
	}
}

string* mapToStrings(std::unordered_map<std::string,bool> Map)
{
	string* names = new string[Map.size()];
	int count =0;
	for ( auto it = Map.begin(); it != Map.end(); ++it, count++ )
    	names[count] = it->first;

 	return names;
}

//generate depth levels of states, from this in all cases  
void generateGlobalQueue(struct state temp, int depth)
{	
	//iterate through the columns & create a new state matrix with the move	

	for (int i=0; i<numcols; i++)
	{
		struct state move = temp;
		if(move.colHeight[i]<6) //means this column has space : its current value is between 0 & 5 inclusive, hence the move in this column is VALID
		{
			int height = move.colHeight[i];
			int pi = abs(depth-3) * 3;
			move.board[height][i] = move.player; //adding the move to the board
			//cout<<"\nPushing the player '"<<move.player<<"'' into column '"<<i<<"' with height '"<<height<<endl;
			move.boardSum+=move.player; //add the player to the boardSum
			move.path[pi++] = height;
			move.path[pi++] = i;
			move.path[pi++] = move.player; //Pushing back this pair into the path taken vector	
			move.player = 3 - move.player; //alternate the player
			move.colHeight[i]++; //add the columnheight

			if(depth == 1)
			{ //if depth = 1, add them to the global queue
				std::string tempHash = giveHash(move);//get The hash of this state 
				//check if this hash already exists in my map
				if (globalHashMap.find(tempHash) != globalHashMap.end()) //which means for find(tempHash) a value was returned that wasn't .end(), hence it exists
					continue;
				else
				{
					globalQueue.push_back(move);
					globalHashMap[tempHash] = 1;
					//globalHashMap.insert({tempHash,1});//add it to the hashmap std::unordered_map<std::string,bool> 
				}
			}		
			else //else call this function back on each of the states
				generateGlobalQueue(move, depth-1);
		}
	}
}

void printQueueOfStates(std::deque<struct state> queueOfStates)
{
	struct state tempQueue;
	int i = 0;

	//iterate through the queue
	for(auto it=queueOfStates.begin(); it!=queueOfStates.end();++it)
    {
		cout<<"\n State "<<++i<<" in the Queue : "<<endl;
		displayBoard(*it);
	}
}

int runAlphaBeta(struct state local, int alpha, int beta, int depth)
{
	struct state temp;
	int tempval = 0;
	temp = local;
	int val;

	if (local.player == 1) //need to maximise, right now its the worst
		val = -9999;
	else
		val = +9999; //need to minimise, right now its the best

	if (depth == 0)
		return evalBoard(local);
	else
	{
		for (int i=0; i<numcols; i++)
		{
			if(local.colHeight[i]<6) //means this column has space : its current value is between 0 & 5 inclusive, hence the move in this column is VALID
			{
				local.board[local.colHeight[i]][i] = local.player; //adding the move to the board
				//local.path.push_back(std::make_pair(std::make_pair (local.colHeight[i],i), local.player)); //Pushing back this pair into the path taken vector	
				local.boardSum+=local.player; //add the player to the boardSum
				local.player = 3 - local.player; //alternate the player
				local.colHeight[i]++; //add the columnheight

				std::string tempHash = giveHash(local);//get The hash of this state 
				//check if this hash already exists in my map
				if (globalHashMap.find(tempHash) != globalHashMap.end()) //which means for find(tempHash) a value was returned that wasn't .end(), hence it exists
					continue;
				else
				{
					globalHashMap[tempHash] = 1;
					tempval = runAlphaBeta(local, alpha, beta, depth -1);
				}

				//tempval = runAlphaBeta(local, alpha, beta, depth -1);

				local = temp;
				if (local.player == 1) //Max player
				{
					if (tempval > val)
						val = tempval; 
				}
				else //Min Player
				{
					if (tempval < val)
						val = tempval;	
				}					
			}

			if (local.player == 1)
			{
				alpha = val;
				if (val > beta) //break if val > Beta
					break;
			}

			else if (local.player == 2)
			{
				beta = val;
				if (val < alpha) //break if val < aplha
					break;
			}	
		}
		return val;
	}
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &pid); //the rank of the current matrix
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs); //total number of processors
	
	int array_of_blocklengths[NBLOCKS] = {12,42,1,1,7};
	MPI_Aint array_of_displacements[NBLOCKS];

	MPI_Datatype array_of_types[NBLOCKS] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
	MPI_Datatype stateType;

	MPI_Aint ex, lowerBound;

	array_of_displacements[0]=0;

	MPI_Type_get_extent(MPI_INT, &lowerBound, &ex);
	array_of_displacements[1] = ex*12;
	array_of_displacements[2] = array_of_displacements[1] + ex*42;
	array_of_displacements[3] = array_of_displacements[2] + ex;
	array_of_displacements[4] = array_of_displacements[3] + ex;

	MPI_Type_struct (5, array_of_blocklengths, array_of_displacements, array_of_types, &stateType);
	MPI_Type_commit (&stateType);

	if(argc < 2) 
	{
        cout<<"You must provide at least one argument\n";
        exit(0);
    }

    //level = argv[1]; //level of the initial state to be generated
    std::string l(argv[1]);
    level = l;

	if (pid == 0)
	{
		if (initialiseState() == -1)
		{
			cout<<"\n Need to Abort";
			exit(0);
		}

		cout<<"\n Processor '"<<pid<<"' : Beginning state - "<<endl;
		displayBoard(startState);

    	cout<<"\n Processor '"<<pid<<"' : Evaluation of the start state : "<<evalBoard(startState)<<endl;
    	
    	generateGlobalQueue(startState, 3); //generating 3 levelled deep states of boards 
		cout<<"\n Processor '"<<pid<<"' : The Global Queue's size : "<<globalQueue.size()<<endl;
		//printQueueOfStates(globalQueue);
	
		globalBestState = startState;
		globalBestVal = -9999;
	}

	MPI_Barrier(MPI_COMM_WORLD); //all threads are here, and the manager has some work ready in its queue

	if (pid == 0)
	{
		starttime = MPI_Wtime();
		cout<<"\n Processor '"<<pid<<"' : Start time :"<<starttime;
		void* sendState;

		int count = globalQueue.size();
		struct state localStartState;
		for (int i=1; i<numProcs; i++) //sending some work for all to start on
		{
			localStartState = globalQueue.at(0);
			globalQueue.pop_front();
			//sendState = void*(localStartState);

			MPI_Send(&localStartState, 1, stateType, i, 0, MPI_COMM_WORLD); //0 corresponds to the state struct
			//MPI_Send(void*(globalHashMap), sizeof(globalHashMap), MPI_DOUBLE, i, 1, MPI_COMM_WORLD); //1 corresponds to the hashmap
			cout<<"\n Sent data for Player "<<localStartState.player<<" with sum "<<localStartState.sum<<" to Proc"<<i;
			
		}

		
		//continuously send the top of the globalQueue until its empty
		while (count)
		{
			int *val;
			//void *receivedState, *receivedHashmap;

			MPI_Status status;
			MPI_Recv (val, 1, MPI_INT, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status); //2 represents the value
			//MPI_Recv (receivedState, sizeof(localStartState), MPI_BYTE, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &status); //3 Represents the state of the matrix
			MPI_Recv(&localStartState, 1, stateType, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &status); //0 corresponds to the state struct
			
			//MPI_RECV (receivedHashmap, 1, MPI_INT, MPI_ANY_SOURCE, 4, MPI_COMM_WORLD, &status); //4 represents the HashMap
			cout<<"\n Received data in 0 for Player "<<localStartState.player<<" with VALUE "<<*val<<" from Proc"<<status.MPI_SOURCE;


			count--;

			//add the hashmap to mine

			localStartState = globalQueue.at(0);

			MPI_Send(&localStartState, 1, stateType, status.MPI_SOURCE, 0, MPI_COMM_WORLD); //0 corresponds to the state struct
			
			//MPI_Send(void*(localStartState), sizeof(localStartState), MPI_BYTE, status.MPI_SOURCE, 0, MPI_COMM_WORLD);//0 corresponds to the state struct
			//MPI_Send(void*(globalHashMap), sizeof(globalHashMap), MPI_DOUBLE, status.MPI_SOURCE, 1, MPI_COMM_WORLD);//1 corresponds to hash map
			
			globalQueue.pop_front();

			if (*val > globalBestVal)
			{
				globalBestState = localStartState;
				globalBestVal = *val;
			}
		}

		//for loop that sends a DUMMY OBJECT with the player value = -1 to all my workers
		for (int i=1; i<numProcs; i++) //sending some work for all to start on
		{
			cout<<"\n Sending an abort message to all procs"<<endl;

			localStartState = startState;
			localStartState.player = -1;

			MPI_Send(&localStartState, 1, stateType, i, 0, MPI_COMM_WORLD); //0 corresponds to the state struct
			//MPI_Send(void*(localStartState), sizeof(localStartState), MPI_BYTE, i, 0, MPI_COMM_WORLD); //0 corresponds to the state struct
			//may not need to send this
				//MPI_ISend(void*(globalHashMap), sizeof(globalHashMap), MPI_DOUBLE, i, 1, MPI_COMM_WORLD); //1 corresponds to the hashmap
		}
	}

	else
	{
		int count = 20;
		cout<<"\n In Processor "<<pid;
		struct state localStartState;
		while (count--) //if I received something from the manager the manager hasnt asked me to quit and ths	
		{
			//void *receivedState, *receivedHashmap;
			cout<<endl<<"Processor "<<pid<<" & count : "<<count;

			MPI_Recv(&localStartState, 1, stateType, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //0 corresponds to the state struct

			cout<<"\n Received data for Player "<<localStartState.player<<" with sum "<<localStartState.boardSum<<" in Proc"<<pid;
			
			//MPI_Recv(receivedState, sizeof(localStartState), MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //0 corresponds to the state struct
			//MPI_Recv(receivedHashmap, sizeof(globalHashMap), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
			//localStartState = struct state *(receivedState);
			if (localStartState.player == -1)
				break;

			//globalHashMap = 
			int value = runAlphaBeta(localStartState, -9999, 9999, 4); 

			MPI_Send (&value, 1, MPI_INT, 0, 2, MPI_COMM_WORLD); //2 represents the value
			MPI_Send(&localStartState, 1, stateType, 0, 3, MPI_COMM_WORLD); //0 corresponds to the state struct
			cout<<"\n Sent data for Player "<<localStartState.player<<" with VALUE "<<value<<" in Proc"<<pid;
			
			
			//MPI_Send (void*(localStartState), sizeof(localStartState), MPI_BYTE,0, 3, MPI_COMM_WORLD); //3 Represents the state of the matrix
			//MPI_Send (void*(globalHashMap), sizeof(globalHashMap), MPI_CHAR, 0, 4, MPI_COMM_WORLD); //4 represents the HashMap
		}
		
	}

	MPI_Barrier(MPI_COMM_WORLD); //post the while loop of the manager, and all the threads are done with their work, hence they're here
	cout<<endl<<"Post While Loop";

	if (pid == 0)
	{
		//call MPI_WTime again
		endtime = MPI_Wtime();
		timeTaken = endtime - starttime;

		//Ultimately print the first pair of the best state, suggesting that path can be taken,
		cout<<"\n The Best first move for you to play would be : ";
		cout<<"Player "<<globalBestState.path[2]<<" : at ["<<globalBestState.path[0]<<","<<globalBestState.path[1]<<"]."<<endl;
		cout<<"\nTime Taken : "<<timeTaken<<endl;

	}

	MPI_Finalize();
 	return 0;
  	
	//displayBoard(globalBestState);
    //displayBoard(startState);
    //testWin();
}


//ADDITIONAL
//If the user wishes to continue from there on, reset the algorithm, with the startState as the current modified State.

/*SINAN's MARSHALL AND DEMARSHALL CODE

void marshall(short *buffer){
buffer[0] = lastCoord.x;
buffer[1] = lastCoord.y;
buffer[2] = lastCoord.z;
buffer[3] = length;
buffer[4] = (short) v;
}

And the demarshall is just a constructor:
Path(short *buffer) : lastCoord(buffer[0],buffer[1],buffer[2]), length(buffer[3]), v((bool) buffer[4]) {}

void marshall(int *buffer)
    {
		int i, j;
		for (i=0; i<path.size(); i++)
			{
				buffer[j++] = path[i].first.first;
				buffer[j++] = path[i].first.second;
				buffer[j++] = path[i].second;
			}
		buffer[j++] = '*';
		buffer[j++] = visited;
		buffer[j++] = player;
		buffer[j++] = boardSum;	

		buffer[j++] = '*';

		//add the matrix board 
		for(int c=0;c<numcols;c++)
			for(int r=0;r<numrows;r++)
				buffer[j++] = board[r][c];

		buffer[j++] = '*';

		//add the column height elements one by one
		for(int c=0;c<colHeight.size();c++)
			buffer[j++] = colHeight[c];

		buffer[j++] = '*';
	}

	state(int *buffer) : lastCoord(buffer[0],buffer[1],buffer[2]), length(buffer[3]), v((bool) buffer[4]) {}

*/





