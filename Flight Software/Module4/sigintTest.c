#include <stdio.h>
#include <unistd.h>
#include <signal.h>

// Global increment value, not sure how to handle this otherwise.
int incr = 1;

// Signal Handler to take any signal
void int_hdl(int signo)
{
    if(signo == SIGINT)
    {
        if(incr < 4)
        {
            incr = incr + 1;
        }
        printf("SIGINT Received! Increment is now %d\n", incr);
    }
    
    if(signo == SIGTERM)
    {
        if(incr > -4)
        {
            incr = incr - 1;
        }
        printf("SIGTERM Received! Increment is now %d\n", incr);
    }
}

int main()
{
    int val = 0;
    int val1 = 0;
    int diff = 0;

    // Compare the particular signals with the handler. If not found, exit.
    if(signal(SIGINT, int_hdl) == SIG_ERR)
    {
        return 1;
    }

    if(signal(SIGTERM, int_hdl) == SIG_ERR)
    {
        return 1;
    }

    while(1)
    {
        val1 = val;
        printf("%d\n", val);
        sleep(1);
        val = val + incr;
        diff = val - val1;
        if (diff >= 2){
            printf("SIGINT Pressed Twice. Exiting Now\n");
            return 1;
        }
    }
    return 0;
}