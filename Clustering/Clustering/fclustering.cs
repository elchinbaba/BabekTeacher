using System;

namespace Clustering
{
	public class FuzzyClustering : Clustering
	{

		public FuzzyClustering(double[,] x) : base(x)
		{
		}

		public FuzzyClustering(double[,] x, int nv) : base(x)
		{
			InitClusters(nv);
		}

		override public void InitClusters()
		{
			this.InitClusters(NC);
		}

		override public void InitClusters(int nv)
		{
			base.InitClusters(nv);
			UpdateMemberships();
		}

		override public void InitClusters(double[][] v)
		{
			base.InitClusters(v);
			UpdateMemberships();
		}

		public double Objective()
		{
			int i, j;
			double s = 0;
			for (i = 0; i < NX; i++)
			{
				for (j = 0; j < NC; j++) s += SqDistanceF(X[i], C[j]) * Math.Pow(U[j][i], M);
			}
			return s;
		}

		// Valid for Euclidean Distance
		public void UpdateMemberships()
		{
			int i, j;
			for (j = 0; j < NC; j++)
			{
				for (i = 0; i < NX; i++)
				{
					double s = 0;
					double d = SqEuclideanDistanceF(X[i], C[j]);
					for (int j1 = 0; j1 < NC; j1++)
					{
						s += Math.Pow(d / SqEuclideanDistanceF(X[i], C[j1]), 1 / (M - 1));
					}
					U[j][i] = 1 / s;
				}
			}
		}

		// Valid for Euclidean Distance
		override public double GetMembership(int c, double[] x)
		{
			double s = 0;
			double d = SqEuclideanDistanceF(x, C[c]);
			if (d == 0) return 1;
			for (int j = 0; j < NC; j++)
			{
				double dj = SqEuclideanDistanceF(x, C[j]);
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
		public void Run(int n)
		{
			for (int it = 0; it < n; it++)
			{
				UpdateMemberships();
				for (int j = 0; j < NC; j++)
				{
					CenterOfGravity(C[j], X, U[j], M);
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
					CenterOfGravity(C[j], X, U[j], M);
				}
				fmin = Objective();
				if (Math.Abs(f - fmin) < eps) break;
			}
			SortClusters();
		}
	}
}
