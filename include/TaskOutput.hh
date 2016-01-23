#ifndef TASKOUTPUT_HH
#define TASKOUTPUT_HH

//MSH_BEGIN
class TaskOutput
{
    public:
        TaskOutput()
        {
            totalDiff=-1;
        }

        TaskOutput(double in1)
        {
            totalDiff=in1;
        }

        TaskOutput& operator=(const TaskOutput& other)
        {
            totalDiff=other.totalDiff;
            return *this;
        }

        double totalDiff; //MSH: primitive
};
//MSH_END
#endif // TASK_OUTPUT_HH

