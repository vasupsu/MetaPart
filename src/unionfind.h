uint64_t findRoot1 (uint64_t *parent, uint64_t elem)
{
	uint64_t root=elem;
	while (parent[root] != root)
	{
		root = parent[root];
//		*numLoads = *numLoads+1;
	}
	while (elem != root)
	{
		uint64_t curr_node = parent[elem];
		parent[elem] = root;
		elem = curr_node;
//		*numLoads = *numLoads+1;
//		*numStores = *numStores+1;
	}
	return elem;
}
uint64_t findRoot (uint8_t *parent, uint64_t elem, int *numRandAccess)
{
	int disp=0;
	uint64_t curParentExtra = (uint64_t)parent[elem*5];
	uint64_t curParent = (uint64_t)(*((uint32_t *)&parent[elem*5 + 1])) + curParentExtra*UINT32_MAX;
	
	if (curParent == 0xBADBADBADl) return elem;
//	if ((elem==14) || (elem==17))
//		disp=1;
        while (curParent != elem)
        {
		//onepass
		if (disp) printf ("%u(%llx)-->%u(%llx)\n", elem, elem, curParent, curParent);
		if (curParent == 0xBADBADBADl) break;
		uint64_t curParentParentExtra = (uint64_t)parent[curParent*5];
		uint64_t curParentParent = (uint64_t)(*((uint32_t *)&parent[curParent*5 +1])) + curParentParentExtra*UINT32_MAX;

		if ((curParent != 0xBADBADBADl) && (curParentParent != 0xBADBADBADl))
		{
			uint32_t tmpParent = (uint32_t)(curParentParent % UINT32_MAX);
			uint8_t tmpRem = (uint8_t) (curParentParent / UINT32_MAX);
			memcpy (&parent[5*elem + 1], &tmpParent, sizeof(uint32_t));
			parent[5*elem] = tmpRem;
			curParent = curParentParent;
//	                parent[elem] = parent[parent[elem]];
		}
		if (curParent != 0xBADBADBADl)
	                elem = curParent;
		else
			break;
		curParentExtra = (uint64_t)parent[elem*5];
		curParent = (uint64_t)(*((uint32_t *)&parent[elem*5 + 1])) + curParentExtra*UINT32_MAX;
                *numRandAccess = *numRandAccess + 1;
//		*numLoads = *numLoads+2;
        }
        return elem;
}

//Assume u<v
void unionOp (uint8_t *parent, uint32_t *sizes, uint64_t u, uint64_t v)
{
	//simple
	uint32_t tmpRem = (uint32_t)(u % UINT32_MAX);
	uint8_t tmpQuot = (uint8_t) u/UINT32_MAX;
	memcpy (&parent[v*5+1], &tmpRem, sizeof(uint32_t));
	parent[v*5] = tmpQuot;
    if (sizes != NULL)
    {
    	sizes[u]+=sizes[v];
	    sizes[v]=0;
    }
//	parent[v]=u;
//	fprintf (fp, "UF %u(size %lu)->%u (size %lu)\n", v, sizes[v], u, sizes[u]);
//	sizes[u]+=sizes[v];

	//weighted
/*	if (sizes[u] < sizes[v])
	{
		parent[u]=v;
		sizes[v]+=sizes[u];
	}
	else
	{
		parent[v]=u;
		sizes[u]+=sizes[v];
	}*/
}
