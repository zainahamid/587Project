/*
Zaina Connect 4 Application to be run in Parallel across Processors

Using MPI & C++

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

7. checkIfVisited() //takes in the state's matrix and returns true if this has been visited. // if it was visited, also makes sure the path's are also copied, is this even required?
8. hash() //to store this particular state in a global hashmap
9. generateNextStep()
10. generateGlobalQueue() : this can choose to have atleast double OR triple the number of 	startStates for the workers to take from the manager 
*/

//#include <mpi.h>
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

int numrows = 6;
int numcols = 7;
double starttime=0;
double endtime=0;
double timeTaken=0;
std::string level;
int numProcs;
std::unordered_map<std::string,bool> globalHashMap; //Maps to the String of the matrix state & the Bool that gives its visited value
std::unordered_map<std::string,bool> localHashMap;

int evaluationTable[6][7] = {{3, 4, 5, 7, 5, 4, 3}, {4, 6, 8, 10, 8, 6, 4}, {5, 8, 11, 13, 11, 8, 5}, {5, 8, 11, 13, 11, 8, 5}, {4, 6, 8, 10, 8, 6, 4}, {3, 4, 5, 7, 5, 4, 3}};
struct state {

    std::vector<std::pair<std::pair<int, int>, int> > path; //vector of pairs 
    bool visited;
    int player;
    double boardSum;
    std::vector<std::vector<int>> board;
    //int **board;
    std::vector<int> colHeight;

    //state() : visited(false), player(1), boardSum(0) {} //structure's constructor
    //write a giveHash() function inside this, that takes the matrix, and sets a hash value for it in the global hash map, 
};

struct state startState;
struct state bestState;
std::deque<struct state> globalQueue;
std::deque<struct state> localQueue;


//global queue of structs state
std::string giveHash(struct state temp)
{
	std::vector<std::vector<int>> board = temp.board;
	std::string hashVal;
    for(int c=0;c<numcols;c++)
    {
		for(int r=0;r<numrows;r++)
		{
			if (board[r][c]==0)
				break;
            else
            	hashVal+=std::to_string(board[r][c]);
		}
		hashVal+="|";
	}		
    return hashVal;
}



void displayBoard(struct state temp)
{
	//displaying the heading of the board
	cout<<"\n\n|";
	for(int col=0; col<numcols; col++)
		cout<<" "<<col+1<<"  ";
	cout<<"|\n";

	std::vector<std::vector<int>> board = temp.board;

	for(int row = numrows-1; row>=0; row--)
	{	cout<<"|";
		for(int col=0; col<numcols; col++)
		{
			if(board[row][col] == 1)
				cout<<" X  ";
			else if (board[row][col] == 2)
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

	cout<<"Column Height : [";
	for(int i =0; i<temp.colHeight.size(); i++)
		cout<<temp.colHeight[i]<<" ";
	cout<<"]"<<endl;;
}


int initialiseState()
{
	startState.player = 1;
	double sum;

	//allocate memory for our state
	std::vector<std::vector<int>> board;

	for(int r=0;r<numrows;r++)
	{
		std::vector<int> row;
		for(int c=0;c<numcols;c++)
			row.push_back(0);

		board.push_back(row);
	}

	if (level.compare("full") == 0)
	{
		board[0][0] = board[0][2] = board[0][3] = board[0][6] = board[1][1] = board[1][2] = board[1][4] = board[1][5]
		= board[2][0] = board[2][4] = board[2][6] = board[3][1] = board[3][5] = board[4][3] = board[4][6] = 1;
		
		board[0][1] = board[0][4] = board[0][5] = board[1][0] = board[1][3] = board[1][6] = board[2][1] = board[2][2]
		= board[2][3] = board[2][5] = board[3][2] = board[3][3] = board[3][4] = board[3][6] = board[4][5] = 2;
		
		startState.boardSum = 45;
	}
	else if (level.compare("middle") == 0)
	{
		board[0][0] = board[0][1] = board[0][2] = board[0][4] = board[1][2] = board[2][3] = board[2][4] = 1;
		board[0][3] = board[1][0] = board[1][1] = board[1][3] = board[1][4] = board[2][2] = board[3][4] = 2;
		startState.boardSum = 21;
	}
	else if (level.compare("begin") == 0)
	{
		board[0][1] = board[0][3] = board[0][4] = 1;
		board[0][2] = board[0][5] = board[0][6] = 2;
		startState.boardSum = 9;
	}
	else //incorrect Input
	{
		cout<<"\n Incorrect Input, Please start again!";
		return -1;
	}

	cout<<"\n In Intialisation"<<endl;

	//set startState.colHeight???
	for(int c=0;c<numcols;c++)
	{
		cout<<"\n Column '"<<c<<"'' : ";
		int height = 0;
		for(int r=0;r<numrows;r++)
		{
			if(board[r][c] != 0)
				cout<<"height = "<<++height<<" & player = "<<board[r][c]<<"."<<endl;
		}
		startState.colHeight.push_back(height);
		cout<<"\n Height of column '"<<c<<"' is '"<<height<<endl;
	}

	startState.board = board;

	cout<<"\n Displaying column of the board 6 : "<<endl;
	for(int r = 0; r<numrows; r++)
		cout<<" "<<board[r][6]<<endl; 

	cout<<"\n Displaying column 6 : "<<endl;
	for(int r = 0; r<numrows; r++)
		cout<<" "<<startState.board[r][6]<<endl; 

	//displayBoard(startState);
	return 0;



	//startState.giveHash(); //visited = giveHash();
}


int checkWin(struct state temp)
{
	std::vector<std::vector<int>> board = temp.board;
	int winner = 0;
	for(int row = 0; row < numrows; row ++)
	{
		for(int col = 0; col < numcols; col ++)
		{
			if(board[row][col] != 0)
			{
				// Across
				if(col < numcols-3)
				{
					if (board[row][col] == board[row][col+1] &&
						board[row][col] == board[row][col+2] &&
						board[row][col] == board[row][col+3])
						{
							winner = board[row][col];
							break;
						}
				}

				// Upwards/downwards
				if(row < numrows-3)
				{
					if (board[row][col] == board[row+1][col] &&
						board[row][col] == board[row+2][col] &&
						board[row][col] == board[row+3][col])
					{
						winner = board[row][col];
						break;
					}
				}
				//Diagonal down-right or up-left
				if(row < numrows-3 && col < numcols-3)
				{
					if (board[row][col] == board[row+1][col+1] &&
						board[row][col] == board[row+2][col+2] &&
						board[row][col] == board[row+3][col+3])
					{
						winner = board[row][col];
						break;
					}
				}
				//Diagonal down-left or up-right
				if(row > 3 && col < numcols-3)
				{
					if (board[row][col] == board[row-1][col+1] &&
						board[row][col] == board[row-2][col+2] &&
						board[row][col] == board[row-3][col+3])
					{
						winner = board[row][col];
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

int testWin ()
{
	std::vector<std::vector<int>> board = startState.board;
	/*board[4][3] = 2;
	board[4][5] = 1;
	board[1][0] = 1;
	board[2][0] = 2;
	board[2][4] = 2;
	board[3][4] = 1;*/
	int winner = checkWin(startState);
	if(winner == 0)
		cout<<"No Winner!"<<endl;
	else
		cout<<"\nWinner = Player "<<winner<<". (1 = X & 2 = O)"<<endl;
	
	displayBoard(startState);
	return 1;
}

int checkDraw(struct state temp)
{
	if((!checkWin(temp)) && temp.boardSum == 63) //when the temp.boardSum == 63 it means the board is filled with 21 1's and 21 2's = 63.
		return 1;
	return 0;
}

int evalBoard(struct state temp)
{
	std::vector<std::vector<int>> board = temp.board;
	int winner = checkWin(temp);
	if (winner == 1)
		return +9999;
	else if (winner == 2)
		return -9999;

	else //No player has won
	{
		if (checkDraw(temp))
			return 0;
		else// AND no draw
		{
			int evalSum = 0;
			for(int r=0;r<numrows;r++)
				for(int c=0;c<numcols;c++)
					if(board[r][c] == 1)
						evalSum += evaluationTable[r][c];
					else if (board[r][c]== 2)
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
			move.board[height][i] = move.player; //adding the move to the board
			cout<<"\nPushing the player '"<<move.player<<"'' into column '"<<i<<"' with height '"<<height<<endl;
			move.boardSum+=move.player; //add the player to the boardSum
			move.player = 3 - move.player; //alternate the player
			
			move.colHeight[i]++; //add the columnheight

			if(depth == 1) //if depth = 1, add them to the global queue
				globalQueue.push_back(move);			
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

int main(int argc, char *argv[])
{
	
	//MPI_Init(&argc, &argv);
	//MPI_Comm_rank(MPI_COMM_WORLD, &rank); //the rank of the current matrix
	//MPI_Comm_size(MPI_COMM_WORLD, &numProcs); //total number of processors
	numProcs = 1;
	if(argc < 2) 
	{
        cout<<"You must provide at least one argument\n";
        exit(0);
    }
    std::string l(argv[1]);
    level = l;

    //level = argv[1]; //level of the initial state to be generated
    if(initialiseState() == -1)
    	exit(0); //fills startState with the type given in level
    displayBoard(startState);
    cout<<"\n Evaluation of the start state : "<<evalBoard(startState)<<endl;


	for(int r = 0;r<numrows;r++)
	{
		cout<<"\n";
		for(int c=0;c<numcols;c++)
		{
			cout<<" "<<startState.board[r][c]<<"  ";
		}
	}

    generateGlobalQueue(startState, 1); //generating 3 levelled deep states of boards 
   
    printQueueOfStates(globalQueue);



    //string* result = mapToStrings(localHashMap); - To conver the local hashmap to a string before sending it to the Manager
    //when the Manager receives this array of strings, it can loop through the array, and add to its global hashmap

    //Have a function that receives HashString value from the worker

    //cout<<"\n Hash Value : "<<giveHash(startState)<<endl;

    /*for (int col = 0; col<numcols; col++)
    {
    	//check by putting startState.player into each of the columns
    	//check if this new state is worth
    		//if yes, add to the queue
    }*/
    
    //from this state, check which player has moves and add all possible moves to a queue of states waiting to be traversed (which don't get pruned)
    //from that queue pick values, and continue adding to the end of the queue those states that can be traversed
    //Each is a queue of structs
    //have a base best state, that keeps getting changed on reaching a terminating state that can be better


    //Ultimately print the first pair of the best state, suggesting that path can be taken, by adding that path to the start state, and displaying start state
    //displayBoard(startState);
   




    
    //testWin();


     //ADDITIONAL
    //If the user wishes to continue from there on, reset the algorithm, with the startState as the current modified State.

}
