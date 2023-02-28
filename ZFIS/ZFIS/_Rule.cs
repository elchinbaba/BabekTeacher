namespace ZFIS
{
	public class _Rule
	{
		internal int[] LeftSideTermIndexes;
		internal int[] RightSideTermIndexes;
		int NLParams;
		int NRParams;

		public _Rule(int[] lterms, int[] rterms)
		{
			//int i;
			NLParams = lterms.Length;
			/*
			LeftSideTermIndexes= new int[NLParams];
			for(i=0; i<NLParams; i++)  LeftSideTermIndexes[i]=lterms[i];
			*/
			NRParams = rterms.Length;
			/*
			RightSideTermIndexes= new int[NRParams];
			for(i=0; i<NRParams; i++)  RightSideTermIndexes[i]=rterms[i];
			*/
			LeftSideTermIndexes = lterms;
			RightSideTermIndexes = rterms;
		}
	}
}
