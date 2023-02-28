namespace DE
{
	public class Optimizer
	{
		internal Problem DEP;

		internal Function CostFun { get { return DEP.CostFun; } }
		internal Function ErrorFun { get { return DEP.ErrorFun; } }

		internal int NPars { get { return DEP.NPars; } }

		internal SearchControls Controls = new SearchControls(60, 0.5, 0.9);

		internal int PopSize { get { return Controls.PopSize; } }
		internal double F { get { return Controls.F; } }
		internal double Cr { get { return Controls.Cr; } }

		internal SearchStrategy Strategy;

		double[] RMin = { 0 };
		double[] RMax = { 1 };

		public Optimizer(Problem dep)
		{
			DEP = dep;
			DEP.DEO = this;
			RNG = new Random();
			Strategy = new BestStrategy(); // default
			Initialized = false;
		}

		double[] Temp; // generated trial
		internal double[] Best; // best ever
		double[] CurBest; // current population's best
		double[][] PNew; // new population
		double[][] PCur; // current population

		internal Random RNG;

		public int Gen = 0; // current generation number

		double MinError;
		double MinCost;

		int BestIndex;

		protected void AllocateMemories()
		{
			Temp = new double[NPars];
			Best = new double[NPars];
			CurBest = new double[NPars];
			PCur = new double[PopSize][];
			PNew = new double[PopSize][];
			for (int i = 0; i < PopSize; i++)
			{
				PCur[i] = new double[NPars];
				PNew[i] = new double[NPars];
			}
		}

		protected bool IsBetter(double c1, double e1, double c2, double e2)
		{
			if (c1 < c2)
			{
				if (e1 <= 0 || e1 <= e2) return true;
				else return false;
			}
			else if (c2 < c1)
			{
				if (e2 <= 0 || e2 <= e1) return false;
				else return true;
			}
			else
			{
				if (e1 < e2) return true;
				else return false;
			}
		}

		protected void InitializePopulation(Ranges ranges)
		{
			AdjustRanges(ranges);
			InitializePopulation();
		}

		protected void InitializePopulation()
		{
			TextOutput.WriteLine("Initializing DE " + NPars + "-D vector space...");
			BestIndex = 0;
			MinCost = MinError = double.MaxValue;
			double[] x;
			double[] trial = DEP.Trial;
			int i = 0;
			if (trial != null)
			{
				x = PCur[0];
				for (int j = 0; j < NPars; j++) x[j] = trial[j];
				MinCost = CostFun(trial);
				MinError = ErrorFun(trial);
				i++;
			}
			for (; i < PopSize; i++)
			{
				x = PCur[i];
				for (int j = 0; j < NPars; j++) x[j] = RNG.Next(RMin[j], RMax[j]);
				double cost = CostFun(x);
				double error = ErrorFun(x);
				if (IsBetter(cost, error, MinCost, MinError))
				{
					MinCost = cost;
					MinError = error;
					BestIndex = i;
				}
			}
			// new: moved from Initialize(...)
			System.Array.Copy(PCur[BestIndex], Best, NPars);
			System.Array.Copy(Best, CurBest, NPars);
			//
			TextOutput.WriteLine("Done.");
		}

		internal void ReInitializePopulation(Ranges ranges)
		{
			if (Initialized) InitializePopulation(ranges);
		}

		internal void ReInitializePopulation()
		{
			if (Initialized) InitializePopulation();
		}

		public bool Initialized { get; protected set; }

		public void Initialize(SearchControls controls, SearchStrategy strategy, Ranges ranges)
		{
			Controls = controls;
			Strategy = strategy;
			AllocateMemories();
			//
			InitializePopulation(ranges);
			Gen = 0;
			strategy.Initialize(RNG);
			Initialized = true;
		}

		/*
		void AdjustRanges(Ranges ranges)
		{
			RMin= ranges.RMin;
			RMax= ranges.RMax;
		}
		*/

		protected internal void AdjustRanges(Ranges ranges)
		{
			double[] rmin = new double[NPars];
			double[] rmax = new double[NPars];
			int n = ranges.RMin.Length;
			if (n <= NPars)
			{
				for (int i = 0; i < n; i++) rmin[i] = ranges.RMin[i];
				for (int i = n; i < NPars; i++) rmin[i] = ranges.RMin[n - 1];
			}
			else
			{
				for (int i = 0; i < NPars; i++) rmin[i] = ranges.RMin[i];
			}
			n = ranges.RMax.Length;
			if (n <= NPars)
			{
				for (int i = 0; i < n; i++) rmax[i] = ranges.RMax[i];
				for (int i = n; i < NPars; i++) rmax[i] = ranges.RMax[n - 1];
			}
			else
			{
				for (int i = 0; i < NPars; i++) rmax[i] = ranges.RMax[i];
			}
			RMin = rmin;
			RMax = rmax;
		}

		protected internal void Run()
		{
			for (int i = 0; i < PopSize; i++)
			{
				System.Array.Copy(PCur[i], Temp, NPars);
				int r1, r2, r3;
				do r1 = RNG.Next(PopSize); while (r1 == i);
				do r2 = RNG.Next(PopSize); while ((r2 == i) || (r2 == r1));
				do r3 = RNG.Next(PopSize); while ((r3 == i) || (r3 == r1) || (r3 == r2));
				//
				Strategy.Apply(F, Cr, NPars, Temp, CurBest, PCur[r1], PCur[r2], PCur[r3]);
				double TestCost = CostFun(Temp);
				double Cost = CostFun(PCur[i]);
				double TestError = ErrorFun(Temp);
				double Error = ErrorFun(PCur[i]);
				if (IsBetter(TestCost, TestError, Cost, Error))
				{
					// substitute i-th solution by test solution
					System.Array.Copy(Temp, PNew[i], NPars);
					if (IsBetter(TestCost, TestError, MinCost, MinError))
					{
						MinCost = TestCost;
						MinError = TestError;
						System.Array.Copy(Temp, Best, NPars);
						BestIndex = i;
						TextOutput.Write("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
						TextOutput.WriteLine((Gen + 1) + "\t" + MinCost + (MinError > 0 ? ("\t" + MinError) : ("")));
					}
				}
				else
				{
					// leave the old solution
					System.Array.Copy(PCur[i], PNew[i], NPars);
				}
			}
			System.Array.Copy(Best, CurBest, NPars);
			// swap population pointers
			double[][] pt = PCur;
			PCur = PNew;
			PNew = pt;
			Gen++;
			TextOutput.Write("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" + Gen);
		}

		public double Run(int iterations, SearchControls controls, SearchStrategy strategy, Ranges ranges)
		{
			Initialize(controls, strategy, ranges);
			return Run(iterations);
		}

		public int SaveFrequency = 0;

		public double Run(int iterations)
		{
			if (!Initialized)
			{
				TextOutput.WriteLine("DE Parameters not set, exited");
				return double.MaxValue;
			}
			TextOutput.WriteLine("Starting DE...");
			TextOutput.WriteLine("Population Size = " + PopSize);
			TextOutput.WriteLine("Initial Cost (and Error) Function value:");
			TextOutput.WriteLine(MinCost + (MinError > 0 ? ("\t" + MinError) : ("")));
			TextOutput.WriteLine("Generation, Cost (and Error) Function Progress:");
			for (; iterations > 0; iterations--)
			{
				if (Controls is AdaptiveControls)
				{
					((AdaptiveControls)Controls).Update();
				}
				Run();
				if (SaveFrequency > 0) if (Gen % SaveFrequency == 0) DEP.Save();
			}
			TextOutput.WriteLine("\nDE stopped.");
			DEP.Save();
			return MinCost;
		}

		public double Run(ExitFunction cond)
		{
			if (!Initialized) return double.MaxValue;
			TextOutput.WriteLine("Starting DE...");
			TextOutput.WriteLine("Population Size = " + PopSize);
			TextOutput.WriteLine("Initial Cost (and Error) Function value:");
			TextOutput.WriteLine(MinCost + (MinError > 0 ? ("\t" + MinError) : ("")));
			TextOutput.WriteLine("Generation, Cost (and Error) Function Progress:");
			while (!cond())
			{
				if (Controls is AdaptiveControls)
				{
					((AdaptiveControls)Controls).Update();
				}
				Run();
				if (SaveFrequency > 0) if (Gen % SaveFrequency == 0) DEP.Save();
			}
			TextOutput.WriteLine("\nDE stopped.");
			DEP.Save();
			return MinCost;
		}
	}
}
