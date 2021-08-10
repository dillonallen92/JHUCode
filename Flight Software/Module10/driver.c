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
    int resourceIDDesired;
    int tickToStart;
    int resourceLockDuration;
    int numTicksLocked;
} Task;


struct Task taskArray[3];
int tick = 0;


int main(int argc, char** argv)
{
    /* Get the command line arguments and place
     * them into the task ids
     */

    if (argc <= 12)
    {
        printf("Not enough input arguments");
    }else{
        /* Command Line Inputs:
         * 1.  Task 1 Priority
         * 2.  Task 1 Resource ID
         * 3.  Task 1 Tick to Start 
         * 4.  Task 1 Tick Count
         * 5.  Task 2 Priority
         * 6.  Task 2 Resource ID
         * 7.  Task 2 Tick to Start
         * 8.  Task 2 Tick Count
         * 9.  Task 3 Priority
         * 10. Task 3 Resource ID
         * 11. Task 3 Tick to Start
         * 12. Task 3 Tick Count
        */
        taskArray[0].originalPriority  = atoi(argv[1]);
        taskArray[0].resourceIDDesired = atoi(argv[2]);
        taskArray[0].tickToStart       = atoi(argv[3]);
        taskArray[0].numTicksLocked    = atoi(argv[4]);
        taskArray[1].originalPriority  = atoi(argv[5]);
        taskArray[1].resourceIDDesired = atoi(argv[6]);
        taskArray[1].tickToStart       = atoi(argv[7]);
        taskArray[1].numTicksLocked    = atoi(argv[8]);
        taskArray[2].originalPriority  = atoi(argv[9]);
        taskArray[2].resourceIDDesired = atoi(argv[10]);
        taskArray[2].tickToStart       = atoi(argv[11]);
        taskArray[2].numTicksLocked    = atoi(argv[12]);


        taskArray[0].currentPriority   = taskArray[0].originalPriority;
        taskArray[1].currentPriority   = taskArray[1].originalPriority;
        taskArray[2].currentPriority   = taskArray[2].originalPriority;


        taskArray[0].resourceIDLocked  = 0;
        taskArray[1].resourceIDLocked  = 0;
        taskArray[2].resourceIDLocked  = 0;

        taskArray[0].resourceLockDuration = taskArray[0].numTicksLocked;
        taskArray[1].resourceLockDuration = taskArray[1].numTicksLocked;
        taskArray[2].resourceLockDuration = taskArray[2].numTicksLocked;

        do
        {
            /*
                * Main Loop Job:
                * 1. call get_next_task_to_run
                * 2. run task on that task
                * 3. update the priority 
                */

            run_task(get_next_task_to_run());
            printf("\nTick %d: Task1 - %d %d %d, Task2 - %d %d %d, Task 3 - %d %d %d\n",
                    tick, taskArray[0].currentPriority, taskArray[0].originalPriority, taskArray[0].resourceIDLocked,
                    taskArray[1].currentPriority, taskArray[1].originalPriority, taskArray[1].resourceIDLocked,
                    taskArray[2].currentPriority, taskArray[2].originalPriority, taskArray[2].resourceIDLocked);        

            tick++;
            update_priority();
        }while(tick <= 40);       
    }
}

void update_priority()
{
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            
            if(taskArray[i].resourceIDLocked == taskArray[j].resourceIDDesired && tick == taskArray[j].tickToStart)
            {
                if(taskArray[i].currentPriority < taskArray[j].currentPriority)
                {
                    taskArray[i].currentPriority = taskArray[j].currentPriority;
                    taskArray[j].tickToStart = taskArray[j].tickToStart + taskArray[i].numTicksLocked;
                }
            }

            if(taskArray[i].resourceIDLocked == 0)
            {
                taskArray[i].currentPriority = taskArray[i].originalPriority;
            }
        }
    }
}

void run_task(int task_id)
{
    int taskArrayID = task_id-1;
    if(taskArray[taskArrayID].resourceIDDesired == 0)
        return;
    else
    {
        taskArray[taskArrayID].resourceIDLocked = taskArray[taskArrayID].resourceIDDesired;
        taskArray[taskArrayID].numTicksLocked--;
        if(taskArray[taskArrayID].numTicksLocked == 0)
        {
            taskArray[taskArrayID].resourceIDLocked  = 0;
            taskArray[taskArrayID].resourceIDDesired = 0;
        } 
    }
}

int get_next_task_to_run()
{
    int task_id = 0;
    for(int i = 0; i < 3; i++)
    {
        if(tick >= taskArray[i].tickToStart && tick <= (taskArray[i].tickToStart + taskArray[i].resourceLockDuration))
        {
            task_id = i+1;
        }
    }
    return task_id;
}
