// @ts-nocheck
// Strict RFC-4180 style parser for EXPORTED scientific reference CSV, not
// a general matrix importer. Unlike the lightweight omics input parser,
// quoted feature IDs may contain commas, CRLF/newlines and escaped quotes.
export function parseReferenceCsv(text){
  if(typeof text!=='string'||!text.length)throw Error('Empty R reference CSV.');
  text=text.replace(/^\uFEFF/,'');
  const rows=[];
  let cells=[],value='',quoted=false,afterQuote=false,started=false;
  function finishCell(){cells.push(value);value='';afterQuote=false;started=false;}
  function finishRow(){finishCell();rows.push(cells);cells=[];}
  for(let i=0;i<text.length;i++){
    const c=text[i];
    if(quoted){
      if(c==='"'&&text[i+1]==='"'){value+='"';i++;}
      else if(c==='"'){quoted=false;afterQuote=true;}
      else value+=c;
      continue;
    }
    if(afterQuote && c!==','&&c!=='\n'&&c!=='\r')
      throw Error('Invalid character after quoted scientific CSV field.');
    if(c==='"'){
      if(started||value.length||afterQuote)
        throw Error('Unescaped scientific CSV quote.');
      quoted=true;started=true;continue;
    }
    if(c===','){finishCell();continue;}
    if(c==='\r'||c==='\n'){
      if(c==='\r'&&text[i+1]==='\n')i++;
      finishRow();continue;
    }
    value+=c;started=true;
  }
  if(quoted)throw Error('Unclosed quoted scientific CSV field.');
  if(started||value.length||afterQuote||cells.length)finishRow();
  if(rows.length<2)throw Error('Reference CSV must have header and fitted rows.');
  const headers=rows.shift();
  if(headers.some(h=>!h)||new Set(headers).size!==headers.length)
    throw Error('Empty or duplicate reference result column names.');
  const data=rows.map((r,i)=>{
    if(r.length!==headers.length)
      throw Error('Ragged scientific CSV row '+(i+2)+'.');
    return Object.fromEntries(headers.map((h,j)=>[h,r[j]]));
  });
  return {headers,rows:data};
}
