using System;

namespace Cluster
{
	public class ProbabilisticClustering : Clustering
	{

		double[][] P; // probabilities of data vectors memberships to clusters P[cluster][data point] - mind the difference with clusterz.cs

		public double[][] Probability { get { return P; } }

		double[] WeightsBuffer; // store vector weights for a cluster

		public ProbabilisticClustering(double[,] x) : base(x)
		{
			WeightsBuffer = new double[NX];
		}

		public ProbabilisticClustering(double[,] x, int nv) : base(x)
		{
			WeightsBuffer = new double[NX];
			InitClusters(nv);
		}

		public double Objective()
		{
			int i, j;
			double s = 0;
			for (i = 0; i < NX; i++)
			{
				for (j = 0; j < NC; j++) s += DistanceF(X[i], C[j]) * Math.Pow(P[j][i], M);
			}
			return s;
		}

		override public void InitClusters()
		{
			this.InitClusters(NC);
		}

		override public void InitClusters(int nv)
		{
			base.InitClusters(nv);
			P = U;
			UpdateMemberships();
		}

		override public void InitClusters(double[][] v)
		{
			base.InitClusters(v);
			P = U;
			UpdateMemberships();
		}

		// Valid for Euclidean Distance
		public void UpdateMemberships()
		{
			int i, j;
			for (i = 0; i < NX; i++)
			{
				for (j = 0; j < NC; j++)
				{
					double p1 = 1;
					for (int k = 0; k < NC; k++)
					{
						if (k != j) p1 *= Math.Pow(EuclideanDistanceF(X[i], C[k]), 1 / (M - 1));
					}
					double s = 0;
					for (int l = 0; l < NC; l++)
					{
						double p2 = 1;
						for (int k = 0; k < NC; k++)
						{
							if (k != l) p2 *= Math.Pow(EuclideanDistanceF(X[i], C[k]), 1 / (M - 1));
						}
						s += p2;
					}
					P[j][i] = p1 / s;
				}
			}
		}

		// Valid for Euclidean Distance
		override public double GetMembership(int c, double[] x)
		{
			double s = 0;
			double d = EuclideanDistanceF(x, C[c]);
			if (d == 0) return 1;
			for (int j = 0; j < NC; j++)
			{
				double dj = EuclideanDistanceF(x, C[j]);
				if (dj == 0) return 0;
				s += Math.Pow(d / dj, 1 / (M - 1));
			}
			return 1 / s;
		}

		double[] x_temp = null;

		override public double GetMembership(int c, int d, double x)
		{
			if (x_temp == null) x_temp = new double[ND];
			for (int k = 0; k < ND; k++)
			{
				if (k == d) x_temp[k] = x;
				else x_temp[k] = Centers[c][k];
			}
			return GetMembership(c, x_temp);
		}

		// Valid for Euclidean Distance
		double[] GetWeights(int c)
		{
			for (int i = 0; i < NX; i++)
				WeightsBuffer[i] = Math.Pow(P[c][i], M) / EuclideanDistanceF(X[i], C[c]);
			return WeightsBuffer;
		}

		// Valid for Euclidean Distance
		public void Run(int n)
		{
			for (int it = 0; it < n; it++)
			{
				UpdateMemberships();
				for (int j = 0; j < NC; j++)
				{
					CenterOfGravity(C[j], X, GetWeights(j));
				}
			}
			SortClusters();
		}

		// Valid for Euclidean Distance
		public void Run(double eps)
		{
			double fmin = Objective();
			while (true)
			{
				double f = fmin;
				UpdateMemberships();
				for (int j = 0; j < NC; j++)
				{
					CenterOfGravity(C[j], X, GetWeights(j));
				}
				fmin = Objective();
				if (Math.Abs(f - fmin) < eps) break;
			}
			SortClusters();
		}

	}
}
