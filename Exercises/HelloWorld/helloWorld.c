#include <stdio.h> // Standard I/O library
#define MAX_USERID_LENGTH 32 // Maximal length of username

// Function for printing "Hello, World!"
void printHelloWorld(void){
	printf("Hello, World!\n");
}

// Gets the username (with max length of MAX_USERID_LENGTH) from the computer/terminal and prints it as "Hello, <USERNAME>!"
void printHelloName(void){
	char username[MAX_USERID_LENGTH];
	cuserid(username);
	printf("Hello, %s!\n", username);	
}

// Main function - Running both exercise A and B
int main(void){
	// Exercise A: Build a hello-world program (a program that outputs "Hello, World!")
	printf("----- Exercise A -----\n");
	printHelloWorld();
	// Exercise B: Make your program also output "Hello, YourUserName!", where YourUserName is your user name (the program must figure it out itself, of course)
	printf("----- Exercise B -----\n");
	printHelloName();
	return 0; // Returning 0 for completed without errors
}
