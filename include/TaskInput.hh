#ifndef TASKINPUT_HH
#define TASKINPUT_HH

//MSH_BEGIN
class TaskInput
{
    public:
        TaskInput()
        {
            g4ndlIndex=-1;
            mcnpIndex=-1;
            sampling=0;
        }

        virtual ~TaskInput()
        {

        }

        TaskInput(int in1, int in2, bool in3)
        {
            g4ndlIndex=in1;
            mcnpIndex=in2;
            sampling=in3;
        }

        TaskInput& operator=(const TaskInput& other)
        {

            g4ndlIndex=other.g4ndlIndex;
            mcnpIndex=other.mcnpIndex;
            sampling=other.sampling;

            return *this;
        }

        int g4ndlIndex; //MSH: primitive
        int mcnpIndex; //MSH: primitive
        bool sampling; //MSH: primitive
};
//MSH_END

#endif // TASK_INPUT_HH
