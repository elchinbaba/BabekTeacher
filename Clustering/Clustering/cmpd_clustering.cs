using System;

namespace Clustering
{
	public class CompoundDEClustering : DE.Problem
	{

		public FuzzyClustering FC { get { return fc; } }
		public ProbabilisticClustering PC { get { return pc; } }

		FuzzyClustering fc;
		ProbabilisticClustering pc;

		int nc, nx, nd;
		//double m;
		double[] xmin, xmax;

		public double[][] Centers { get { return fc.Centers; } }

		public int NCusters { get { return nc; } }
		public int NData { get { return nx; } }
		public int Dim { get { return nd; } }

		public double PFuzzifier { get { return pc.Fuzzifier; } set { pc.Fuzzifier = value; pc.UpdateMemberships(); } }
		public double FFuzzifier { get { return fc.Fuzzifier; } set { fc.Fuzzifier = value; fc.UpdateMemberships(); } }

		public CompoundDEClustering(double[,] x, int nc)
		{
			fc = new FuzzyClustering(x);
			pc = new ProbabilisticClustering(x);
			this.nc = nc;
			fc.InitClusters(nc);
			pc.InitClusters(nc);
			nd = fc.Dim;
			nx = x.Length / nd;
			//m= pc.Fuzzifier;
			xmin = new double[nc * nd];
			xmax = new double[nc * nd];
			Trial = new double[nc * nd];
			int l = 0;
			for (int j = 0; j < nc; j++) for (int k = 0; k < nd; k++) this.Trial[l++] = fc.Centers[j][k];
			l = 0;
			for (int j = 0; j < nc; j++) for (int k = 0; k < nd; k++) { xmin[l] = fc.XMinBounds[k]; xmax[l] = fc.XMaxBounds[k]; l++; }
			NPars = nc * nd; // DE.Problem dimension
		}

		public void InitClusters(double[][] v)
		{
			fc.InitClusters(v);
			pc.InitClusters(v);
			this.nc = fc.NClusters;
			//
			int l = 0;
			for (int j = 0; j < nc; j++) for (int k = 0; k < nd; k++) this.Trial[l++] = v[j][k];
			//
			fc.UpdateMemberships();
			pc.UpdateMemberships();
		}

		public double Objective()
		{
			return fc.Objective() + pc.Objective();
		}

		override public double EvaluateCost(double[] pars)
		{
			int l = 0;
			for (int j = 0; j < nc; j++) for (int k = 0; k < nd; k++) fc.Centers[j][k] = pc.Centers[j][k] = pars[l++];
			fc.UpdateMemberships();
			pc.UpdateMemberships();
			return Objective();
		}

		bool de_active_ = false;
		DE.SearchControls de_controls_;
		DE.SearchStrategy de_strategy_;

		void Run_(int steps)
		{
			Optimize(steps, de_controls_, de_strategy_, new DE.Ranges(xmin, xmax));
			int l = 0;
			for (int j = 0; j < nc; j++) for (int k = 0; k < nd; k++) fc.Centers[j][k] = pc.Centers[j][k] = Solution[l++];
			//
			fc.UpdateMemberships();
			pc.UpdateMemberships();
			SortClusters();
		}

		public void Run(int steps, int ps)
		{
			de_active_ = true;
			de_controls_ = new DE.SearchControls(ps, 0.9, 1.0);
			de_strategy_ = new DE.RandStrategy();
			Run_(steps);
		}

		public void SortClusters()
		{
			fc.SortClusters();
			pc.SortClusters();
		}

		public void Run(int steps)
		{
			if (!de_active_) Run(steps, 100);
			else
			{
				Run_(steps);
			}
		}

		public void OutputClusterMDegree(int c)
		{
			for (int i = 0; i < nx; i++) Console.WriteLine(fc.Membership[c][i]); // membership degree
		}

		public void OutputClusterMProbability(int c)
		{
			for (int i = 0; i < nx; i++) Console.WriteLine(pc.Probability[c][i]); // membership probability
		}

		public double[][] MembershipDegree { get { return fc.Membership; } }
		public double[][] MembershipProbability { get { return pc.Probability; } }

	}
}
