#pragma once
// from http://roar11.com/2016/01/a-platform-independent-thread-pool-using-c14/

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <utility>
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <functional>
#include <future>
#include <memory>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace Threading
{
	template <typename T>
	class ThreadSafeQueue
	{
	public:

		~ThreadSafeQueue(void)
		{
			invalidate();
		}

		/**
		* Attempt to get the first value in the queue.
		* Returns true if a value was successfully written to the out parameter, false otherwise.
		*/
		bool tryPop(T& out)
		{
			std::lock_guard<std::mutex> lock{ m_mutex };
			if (m_queue.empty() || !m_valid)
			{
				return false;
			}
			out = std::move(m_queue.front());
			m_queue.pop();
			return true;
		}

		/**
		* Get the first value in the queue.
		* Will block until a value is available unless clear is called or the instance is destructed.
		* Returns true if a value was successfully written to the out parameter, false otherwise.
		*/
		bool waitPop(T& out)
		{
			std::unique_lock<std::mutex> lock{ m_mutex };
			m_condition.wait(lock, [this]()
			{
				return !m_queue.empty() || !m_valid;
			});
			/*
			* Using the condition in the predicate ensures that spurious wakeups with a valid
			* but empty queue will not proceed, so only need to check for validity before proceeding.
			*/
			if (!m_valid)
			{
				return false;
			}
			out = std::move(m_queue.front());
			m_queue.pop();
			return true;
		}

		/**
		* Push a new value onto the queue.
		*/
		void push(T value)
		{
			std::lock_guard<std::mutex> lock{ m_mutex };
			m_queue.push(std::move(value));
			m_condition.notify_one();
		}

		/**
		* Check whether or not the queue is empty.
		*/
		bool empty(void) const
		{
			std::lock_guard<std::mutex> lock{ m_mutex };
			return m_queue.empty();
		}

		/**
		* Clear all items from the queue.
		*/
		void clear(void)
		{
			std::lock_guard<std::mutex> lock{ m_mutex };
			while (!m_queue.empty())
			{
				m_queue.pop();
			}
			m_condition.notify_all();
		}

		/**
		* Invalidate the queue.
		* Used to ensure no conditions are being waited on in waitPop when
		* a thread or the application is trying to exit.
		* The queue is invalid after calling this method and it is an error
		* to continue using a queue after this method has been called.
		*/
		void invalidate(void)
		{
			std::lock_guard<std::mutex> lock{ m_mutex };
			m_valid = false;
			m_condition.notify_all();
		}

		/**
		* Returns whether or not this queue is valid.
		*/
		bool isValid(void) const
		{
			std::lock_guard<std::mutex> lock{ m_mutex };
			return m_valid;
		}

	private:
		std::atomic_bool m_valid{ true };
		mutable std::mutex m_mutex;
		std::queue<T> m_queue;
		std::condition_variable m_condition;
	};

	class ThreadPool
	{
	private:
		class IThreadTask
		{
		public:
			IThreadTask(void) = default;
			virtual ~IThreadTask(void) = default;
			IThreadTask(const IThreadTask& rhs) = delete;
			IThreadTask& operator=(const IThreadTask& rhs) = delete;
			IThreadTask(IThreadTask&& other) = default;
			IThreadTask& operator=(IThreadTask&& other) = default;

			virtual void execute() = 0;
		};

		template <typename Func>
		class ThreadTask : public IThreadTask
		{
		public:
			ThreadTask(Func&& func) : m_func{ std::move(func) }
			{
			}

			~ThreadTask(void) override = default;
			ThreadTask(const ThreadTask& rhs) = delete;
			ThreadTask& operator=(const ThreadTask& rhs) = delete;
			ThreadTask(ThreadTask&& other) = default;
			ThreadTask& operator=(ThreadTask&& other) = default;

			void execute() override
			{
				m_func();
			}

		private:
			Func m_func;
		};

	public:
		/**
		* A wrapper around a std::future that adds the behavior of futures returned from std::async.
		* Specifically, this object will block and wait for execution to finish before going out of scope.
		*/
		template <typename T>
		class TaskFuture
		{
		public:
			TaskFuture(std::future<T>&& future)
				:m_future{ std::move(future) }
			{
			}

			TaskFuture(const TaskFuture& rhs) = delete;
			TaskFuture& operator=(const TaskFuture& rhs) = delete;
			TaskFuture(TaskFuture&& other) = default;
			TaskFuture& operator=(TaskFuture&& other) = default;
			~TaskFuture(void)
			{
				if (m_future.valid())
				{
					m_future.get();
				}
			}

			auto get(void)
			{
				return m_future.get();
			}


		private:
			std::future<T> m_future;
		};

	public:

		explicit ThreadPool(int numThreads = -1)
			:m_done{ false },
			m_workQueue{},
			m_threads{}
		{
			numThreads = numThreads > 0 ? numThreads : std::max(std::thread::hardware_concurrency() - 1u, 1u);
			try
			{
				for (int i = 0; i < numThreads; ++i)
				{
					m_threads.emplace_back(&ThreadPool::worker, this);
				}
			}
			catch (...)
			{
				destroy();
				throw;
			}
		}
		~ThreadPool(void)
		{
			destroy();
		}

		ThreadPool(const ThreadPool& rhs) = delete;
		ThreadPool& operator=(const ThreadPool& rhs) = delete;
	
		template <typename Func, typename... Args>
		auto submit(Func&& func, Args&&... args)
		{
			auto boundTask = std::bind(std::forward<Func>(func), std::forward<Args>(args)...);
			using ResultType = std::invoke_result_t<std::decay_t<decltype(boundTask)>>;
			//using ResultType = std::result_of_t<decltype(boundTask)()>;
			using PackagedTask = std::packaged_task<ResultType()>;
			using TaskType = ThreadTask<PackagedTask>;

			PackagedTask task{ std::move(boundTask) };
			TaskFuture<ResultType> result{ task.get_future() };
			m_workQueue.push(std::make_unique<TaskType>(std::move(task)));
			return result;
		}

	private:
		void worker(void)
		{
			// Constantly running function each thread uses to acquire work items from the queue.
			while (!m_done)
			{
				std::unique_ptr<IThreadTask> pTask{ nullptr };
				if (m_workQueue.waitPop(pTask))
				{
					pTask->execute();
				}
			}
		}

		void destroy(void)
		{
			// Invalidates the queue and joins all running threads.
			m_done = true;
			m_workQueue.invalidate();
			for (auto& thread : m_threads)
			{
				if (thread.joinable())
				{
					thread.join();
				}
			}
		}

	private:
		std::atomic_bool m_done;
		ThreadSafeQueue<std::unique_ptr<IThreadTask>> m_workQueue;
		std::vector<std::thread> m_threads;
	};

	//namespace DefaultThreadPool
	//{
	//	/**
	//	* Get the default thread pool for the application.
	//	* This pool is created with std::thread::hardware_concurrency() - 1 threads.
	//	*/
	//	inline ThreadPool& getThreadPool(void)
	//	{
	//		static ThreadPool defaultPool;
	//		return defaultPool;
	//	}

	//	/**
	//	* Submit a job to the default thread pool.
	//	*/
	//	template <typename Func, typename... Args>
	//	inline auto submitJob(Func&& func, Args&&... args)
	//	{
	//		return getThreadPool().submit(std::forward<Func>(func), std::forward<Args>(args)...);
	//	}
	//}

}

// msmpi.lib;
