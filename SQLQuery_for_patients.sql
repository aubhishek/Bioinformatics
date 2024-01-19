select  C1.PatientDurableKey, C1.ProcedureTerminologyCode, C1.ProcedureTerminologyName, C1.ProcedureTerminologyCodeSet from CDW_NEW.deid_uf.CodedProcedureFact C1
where C1.PatientDurableKey in (select PatientDurableKey
from (
    select PatientDurableKey, value, referencevalues,
	       case when upper_limit is null then 2.8
		        else upper_limit
	       end as upper_limit_value
    from (
	    select PatientDurableKey, value, referencevalues, 
	           stuff(stuff(ref_str+'x', patindex('%[0-9][^0-9.]%', ref_str+'x') + 1, len(ref_str), ''
                  ), 1, patindex('%[0-9]%', ref_str) - 1, '') as upper_limit
	    from (
	        select l2.PatientDurableKey, l2.value, l2.referencevalues,  
	               right(l2.referencevalues, charindex(':', reverse(l2.referencevalues) + ':') - 1) as ref_str
            from CDW_NEW.deid_uf.LabComponentResultFact l2
            where  l2.ProcedureName like 'Porphobil%urine%'
            and l2.ComponentName like '%porphobilinogen%urine%'
            and l2.NumericValue is not null 
	    ) tmp_1
    ) tmp_2
) tmp_3
where cast (value as float) < 2 * cast (upper_limit_value as float) 
union
select distinct PatientDurableKey
from ( 
    select PatientDurableKey, value, referencevalues,
	       case when upper_limit is null then 6.2
		        else upper_limit
	       end as upper_limit_value
    from (
	    select PatientDurableKey, value, referencevalues, 
	           stuff(stuff(ref_str+'x', patindex('%[0-9][^0-9.]%', ref_str+'x') + 1, len(ref_str), ''
                  ), 1, patindex('%[0-9]%', ref_str) - 1, '') as upper_limit
	    from (
	        select l2.PatientDurableKey, l2.value, l2.referencevalues,  
	               right(l2.referencevalues, charindex(':', reverse(l2.referencevalues) + ':') - 1) as ref_str
            from CDW_NEW.deid_uf.LabComponentResultFact l2
            where  l2.ProcedureName like '%levulinic%urine%'
            and l2.ComponentName like '%aminolevulinic%urine%'
            and l2.NumericValue is not null 
	    ) tmp_1
    ) tmp_2
) tmp_3
where cast (value as float) < 2 * cast (upper_limit_value as float) )