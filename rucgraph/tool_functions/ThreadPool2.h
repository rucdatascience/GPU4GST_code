/*this is from https://github.com/progschj/ThreadPool2 

an explanation: https://www.cnblogs.com/chenleideblog/p/12915534.html

ThreadPool uses mutex, while ThreadPool2 uses shared_mutex

std::shared_mutex��std::mutex�����ܶԱ�(benchmark):
https://blog.csdn.net/analogous_love/article/details/97918304
*/


#ifndef THREAD_POOL_H2
#define THREAD_POOL_H2

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <shared_mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

class ThreadPool2 {
public:
    ThreadPool2(size_t); // �̳߳صĹ��캯��
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) 
        -> std::future<typename std::result_of<F(Args...)>::type>; // ��������ӵ��̳߳ص����������
    ~ThreadPool2(); // �̳߳ص���������
private:
    // need to keep track of threads so we can join them
    std::vector< std::thread > workers; // ���ڴ���̵߳����飬��vector��������
    // the task queue
    std::queue< std::function<void()> > tasks; // ���ڴ������Ķ��У���queue���н��б��档��������Ϊstd::function<void()>����Ϊstd::function��ͨ�ö�̬������װ������������������д�ŵ���һ��������
    
    // synchronization
    std::shared_mutex queue_mutex; // һ������������еĻ��������ڲ�����������߳�ȡ��������Ҫ�������������а�ȫ����

    /*�˴�������condition_variable����Ϊ����רΪmutex��Ƶģ���������shared_mutex*/
    std::condition_variable_any condition; // һ������֪ͨ�߳��������״̬����������������������֪ͨ�߳̿���ִ�У��������wait״̬

    bool stop; // ��ʶ�̳߳ص�״̬�����ڹ����������ж��̳߳�״̬���˽�
};
 
// the constructor just launches some amount of workers
inline ThreadPool2::ThreadPool2(size_t threads) // ���캯������Ϊinline�����ղ���threads��ʾ�̳߳���Ҫ�������ٸ��̡߳�
    :   stop(false) // stop��ʼΪfalse������ʾ�̳߳������š�
{
    for(size_t i = 0;i<threads;++i) // ����forѭ�������δ���threads���̣߳��������߳�����workers�С�
        workers.emplace_back( 
            /*
            ��vector�У�emplace_back()��Ա������������������β������һ����������Ч����push_back()һ����������������΢���죬
            ��emplace_back(args)�з���Ķ���Ĳ�������push_back(OBJ(args))�з�����Ƕ��󡣼�emplace_back()ֱ�����������Դ���Ĳ���ֱ�ӵ��ö���Ĺ��캯�������µĶ���
            ��push_back()���ȵ��ö���Ĺ��캯������һ����ʱ�����ٽ���ʱ���󿽱��������ڴ��С�
            
            lambda���ʽ�ĸ�ʽΪ��
            [ ���� ] ( �β� ) ˵����(��ѡ) �쳣˵�� attr -> �������� { ������ }
            ��������lambda���ʽΪ [ ���� ] { ������ } ���͡������lambda���ʽ���Ǵ����̵߳�ִ�к���.       
            */
            [this] // ��lambda���ʽ�����̳߳�ָ��this�����ں�������ʹ�ã������̳߳س�Ա����stop��tasks�ȣ�    
            {
                for(;;) // for(;;)Ϊһ����ѭ������ʾÿ���̶߳��ᷴ������ִ�У�����ʵÿ���̳߳��е��̶߳���������
                {
                    std::function<void()> task; // ��ѭ���У��ȴ���һ����װvoid()������std::function����task�����ڽ��պ�������������е�������ʵ����

                    {
                        std::unique_lock<std::shared_mutex> lock(this->queue_mutex); // ��{}�ڣ�queue_mutex������״̬

                        /*�����ʾ���̳߳���ֹͣ������������в�Ϊ�գ��򲻻���뵽wait״̬��
                          ���ڸտ�ʼ�����̳߳أ��̳߳ر�ʾδֹͣ�����������Ϊ�գ�����ÿ���̶߳�����뵽wait״̬��
                          ������������������֪ͨ���߳̾ͻ�������½���
                        */
                        this->condition.wait(lock,
                            [this]{ return this->stop || !this->tasks.empty(); });

                        /*���̳߳��Ѿ�ֹͣ���������Ϊ�գ���return���������߳�������ѭ�������׳��̳߳أ��ڴ�֮ǰ��ÿ���߳�Ҫô��wait״̬��Ҫô��ִ�������task*/
                        if(this->stop && this->tasks.empty())
                            return;

                        /*
                        ����������еĵ�һ��������task��ǣ�Ȼ����������и����񵯳������˴��߳�ʵ�ڻ������������еĻ�����������½��еģ����������������̺߳�
                        �߳���wait�����ڵõ���������еĻ������Ż��������ִ�С���������ֻ����һ���߳��õ����񣬲��ᷢ����ȺЧӦ��
                        ���˳���{ }�����Ƕ�������е����ӵ���Ҳ�ͷ��ˣ�Ȼ�����ǵ��߳̾Ϳ���ִ�������õ�������task�ˣ�ִ�����֮���߳��ֽ�������ѭ����
                        */
                        task = std::move(this->tasks.front());
                        this->tasks.pop(); // task����ѻ���queue��
                    }

                    task();
                }
            }
        );
}

// add new work item to the pool
/*
equeue��һ��ģ�庯�����������β�ΪF��Args������class... Args��ʾ��������βΡ�
auto�����Զ��Ƶ���equeue�ķ������ͣ��������β�Ϊ(F&& f, Args&&... args)������&&��ʾ��ֵ���á���ʾ����һ��F���͵�f�������ɸ�Args���͵�args��

typename std::result_of<F(Args...)>::type   //�����ArgsΪ������F�ĺ������͵ķ�������
std::future<typename std::result_of<F(Args...)>::type> //std::future���������첽�����Ľ��
���շ��ص��Ƿ���std::future�е�F(Args��)�������͵��첽ִ�н����
*/
template<class F, class... Args>
auto ThreadPool2::enqueue(F&& f, Args&&... args) 
    -> std::future<typename std::result_of<F(Args...)>::type> // ��ʾ�������ͣ���lambda���ʽ�еı�ʾ����һ����
{
    using return_type = typename std::result_of<F(Args...)>::type; // �����ArgsΪ������F�ĺ������͵ķ�������

    auto task = std::make_shared< std::packaged_task<return_type()> >(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
        
    std::future<return_type> res = task->get_future(); // res�б���������Ϊreturn_type�ı�������task�첽ִ����ϲſ��Խ�ֵ�����ȥ
    {
        std::unique_lock<std::shared_mutex> lock(queue_mutex);

        // don't allow enqueueing after stopping the pool
        if(stop)
            throw std::runtime_error("enqueue on stopped ThreadPool2");

        tasks.emplace([task](){ (*task)(); });
    }
    condition.notify_one(); // //�������������к���Ҫȥ����һ���߳�
    return res; // //���߳�ִ����ϣ���ִ�еĽ������
}

// the destructor joins all threads
inline ThreadPool2::~ThreadPool2()
{
    /*
    �����������У��ȶ���������м�������ֹͣ�������Ϊtrue������������ʹ���µĲ����������Ҳ��ִ��ʧ��
    */
    {
        std::unique_lock<std::shared_mutex> lock(queue_mutex); // ��{}�ڣ�queue_mutex������״̬
        stop = true;
    }
    /*ʹ�������������������̣߳������̶߳�������ִ��*/
    condition.notify_all();

    /*��stop����Ϊtrue�����������Ϊ��ʱ����Ӧ���߳̽�������ѭ������
    ��ÿ���߳�����Ϊjoin���ȵ�ÿ���߳̽�����Ϻ����߳����˳���*/
    for(std::thread &worker: workers)
        worker.join();
}

#endif











/*Example:

it seems that ThreadPool is faster than ThreadPool2 on Linux: 0.002s vs 0.003s,
but is slower than ThreadPool2 on Windows: 0.005s vs 0.002s,
for num = 1e3 in the following codes.

-------------------------------

#include <iostream>
#include <tool_functions/ThreadPool.h>
#include <tool_functions/ThreadPool2.h>
using namespace std;

double func(int i) {
    return i;
}

void ThreadPool1_2_compare()
{
    int num = 1e3;
    {
        auto begin = std::chrono::high_resolution_clock::now();
        ThreadPool pool(4); // ����һ���̳߳أ������߳�Ϊ4
        std::vector<std::future<double>> results;
        for (int i = 0; i < num; ++i)
        {
            results.emplace_back(
                pool.enqueue([i] {
                    return func(i);
                    }));
        }
        for (auto&& result : results)
            result.get();
        auto end = std::chrono::high_resolution_clock::now();
        double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
        cout << "ThreadPool: " << runningtime << "s" << endl;
    }

    {
        auto begin = std::chrono::high_resolution_clock::now();
        ThreadPool2 pool(4); // ����һ���̳߳أ������߳�Ϊ4
        std::vector<std::future<double>> results;
        for (int i = 0; i < num; ++i)
        {
            results.emplace_back(
                pool.enqueue([i] {
                    return func(i);
                    }));
        }
        for (auto&& result : results)
            result.get();
        auto end = std::chrono::high_resolution_clock::now();
        double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
        cout << "ThreadPool2: " << runningtime << "s" << endl;
    }



}

int main()
{
    ThreadPool1_2_compare();
}

*/



