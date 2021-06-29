#include <stdio.h>
#include <unistd.h>
#include <signal.h>

void int_hdl(int signo)
{
    if(signo == SIGINT)
    {
        printf("pressed the SIGINT thing\n");
    }
    
    if(signo == SIGTERM)
    {
        printf("Pressed the SIGTERM thing\n");
    }
}

int main()
{
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
        printf("Hello world\n");
        sleep(1);
    }
    return 0;
}