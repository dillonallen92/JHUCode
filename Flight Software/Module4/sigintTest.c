#include <stdio.h>
#include <unistd.h>
#include <signal.h>

// Global increment value, not sure how to handle this otherwise.
int val = 0;
int incr = 1;

// Global flag variable? Trying a hint from teammate
int sigintFlag = 0;

// code snippet from stack overflow, adopted to my use
// link: https://stackoverflow.com/questions/21542077/c-sigalrm-alarm-to-display-message-every-second
// Varied the function by doing the arithmetic and printing the values
void display_message(int s)
{
    printf("\n%d\n", val);
    val = val + incr;
    sigintFlag = 0;
    alarm(1);
    signal(SIGALRM, display_message);
}

// Signal Handler to take any signal
void int_hdl(int signo)
{
    if(signo == SIGINT)
    {
        sigintFlag = sigintFlag + 1;
        if (sigintFlag > 1)
        {
            printf("\nExiting program now.\n");
        }else
        {
            if(incr < 4)
            {
                incr = incr + 1;
            }
            printf("\nSIGINT Received! Increment is now %d\n", incr);
        }
        
    }
    
    if(signo == SIGTERM)
    {
        if(incr > -4)
        {
            incr = incr - 1;
        }
        printf("\nSIGTERM Received! Increment is now %d\n", incr);
    }
}

int main()
{
    int val = 0;

    // Compare the particular signals with the handler. If not found, exit.
    if(signal(SIGINT, int_hdl) == SIG_ERR)
    {
        return 1;
    }

    if(signal(SIGTERM, int_hdl) == SIG_ERR)
    {
        return 1;
    }

    // noticed stuff would only print if I had a "\n" preceeding the print in my handler
    // I wondered if a flush would work to start printing, hence this is here. I am pretty sure this
    // isnt proper code.
    //fflush(stdout);

    // The only other way I found on stack overflow without using sleep(n)
    signal(SIGALRM, display_message);
    alarm(1);
    while(1)
    {
        if (sigintFlag > 1)
        {
            return 1;
        }
    }
    return 0;
}