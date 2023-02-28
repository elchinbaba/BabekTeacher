using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DE
{
	public class Problem
	{
		protected internal double[] Trial = null;

		public double[] Solution { get { return DEO.Best; } }

		protected internal int NPars = 1;

		protected internal Function CostFun;
		protected internal Function ErrorFun;

		public Problem(Function costf, Function errorf, double[] trial) : this()
		{
			CostFun = costf;
			ErrorFun = errorf;
			if (trial != null)
			{
				Trial = trial;
				NPars = Trial.Length;
			}
			else
			{
				NPars = 1;
				Trial = null;
			}
		}

		public Problem(Function costf, double[] trial) : this(costf, delegate (double[] x) { return 0; }, trial)
		{ }

		public Problem(Function costf, Function errorf, int npars) : this()
		{
			CostFun = costf;
			ErrorFun = errorf;
			Trial = null;
			NPars = npars;
			if (NPars < 0) NPars = 1;
		}

		public Problem(Function costf, int npars) : this(costf, delegate (double[] x) { return 0; }, npars)
		{ }

		protected Problem()
		{
			CostFun = EvaluateCost;
			ErrorFun = EvaluateError;
			DEO = new Optimizer(this);
			// needs to be initialized in derived class' constructor
			// override EvaluateCost(.) (and EvaluateError(.) and Save(.), if necessary)
			// set value of NPars
		}

		virtual public double EvaluateCost(double[] temp) { return 0; }

		virtual public double EvaluateError(double[] temp) { return 0; }

		virtual public int Save() { return 1; }

		public Optimizer DEO;

		public void Initialize(SearchControls controls, SearchStrategy strategy, Ranges ranges)
		{
			DEO.Initialize(controls, strategy, ranges);
		}

		public double Optimize(int iterations)
		{
			return DEO.Run(iterations);
		}

		public double Optimize(int iterations, SearchControls controls, SearchStrategy strategy, Ranges ranges)
		{
			DEO.Initialize(controls, strategy, ranges);
			return DEO.Run(iterations);
		}

		public double Optimize(ExitFunction cond)
		{
			return DEO.Run(cond);
		}

		public bool IsFeasible { get { return Solution == null || ErrorFun(Solution) > 0.0 ? false : true; } }
		public double Cost { get { return Solution == null ? double.MaxValue : CostFun(Solution); } }
		public double Error { get { return Solution == null ? double.MaxValue : ErrorFun(Solution); } }

		public void ResetSolution()
		{
			DEO.ReInitializePopulation();
		}

		public void ResetSolution(Ranges ranges)
		{
			DEO.ReInitializePopulation(ranges);
		}

		public void ResetSolution(double[] trial)
		{
			Trial = trial;
			DEO.ReInitializePopulation();
		}

		public void ResetSolution(Ranges ranges, double[] trial)
		{
			Trial = trial;
			DEO.ReInitializePopulation(ranges);
		}
	}
}
