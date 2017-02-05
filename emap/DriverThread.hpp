#pragma once

#include "readerwriterqueue/readerwriterqueue.h"

#include <thread>

// convenience class, provides queueing between threads for speed
// Data objects are moved between worker and device threads,
// so worker may send commands back to device too
template <typename Data>
class DriverThread
{
protected:
	// return false when last data object is filled (device thread)
	virtual bool fill(Data & data) = 0;

	// return false when last data object is processed (worker thread)
	virtual bool process(Data & data) = 0;

	template<typename F> DriverThread(F newData)
	: outgoingQueue(16), reuseQueue(16)
	{
		for (int i = 0; i < 16; ++ i)
			reuseQueue.enqueue(std::move(newData()));
	}

	// call in constructor
	void start()
	{
		workTh = std::thread(&DriverThread<Data>::workRunner, this);
		devTh = std::thread(&DriverThread<Data>::devRunner, this);
	}

	// call in destructor
	void join()
	{
		devTh.join();
		workTh.join();
	}

private:
	std::thread workTh;
	std::thread devTh;

	moodycamel::BlockingReaderWriterQueue<Data> outgoingQueue;
	moodycamel::BlockingReaderWriterQueue<Data> reuseQueue;

	void devRunner()
	{	
		Data data;
		for (;;) {
			reuseQueue.wait_dequeue(data);
			if (fill(data)) {
				// TODO: at least spew a message out when buffers are dropped
				outgoingQueue.try_enqueue(std::move(data));
			} else {
				outgoingQueue.enqueue(std::move(data));
				break;
			}
		}
	}

	void workRunner()
	{
		Data data;
		for (;;) {
			outgoingQueue.wait_dequeue(data);
			if (!process(data))
				break;
			reuseQueue.enqueue(std::move(data));
		}
	}
};
