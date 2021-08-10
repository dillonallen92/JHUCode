#include <stdio.h>

/*
* Task Management Program
*  Implement a simple priority-inheritance algorithm and test it in a sample application. Your program will manage
*  3 tasks (T1, T2, and T3) which can lock one of three resources (R1, R2, R3) for a period of time measured in ticks.
*  The system implements 3 priorities (0 = low, 1 = medium, 2 = high).
*
*/

void update_priority();
void run_task(int task_id);
int  get_next_task_to_run();

typedef struct Task
{
    int taskID;
    int originalPriority;
    int currentPriority;
    int resourceIDLocked;
    int tick;
    int resourceLockedTicks;
} Task;

struct Task taskArray[3];


int main()
{
    taskArray[0].taskID = 1;
    taskArray[1].taskID = 2;
    taskArray[2].taskID = 3;
    update_priority();
}

void update_priority()
{
    for(int i = 0; i < 3; i++)
    {
        printf("Task ID: %d\n",taskArray[i].taskID);
    }
}

void run_task(int task_id)
{

}