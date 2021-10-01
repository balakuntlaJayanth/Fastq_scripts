package main

import (
  "bufio"
  "bytes"
  "fmt"
  "io"
  "os"
)

func main() {
  n, sLen, qLen := 0, int64(0), int64(0)
  q10_20 := int64(0)
  q20_30 := int64(0)
  q30 := int64(0)
  q10 := int64(0)
  var fqr FqReader
  fqr.Reader = bufio.NewReader(os.Stdin)
  for r, done := fqr.Iter(); !done; r, done = fqr.Iter() {
    n += 1
    sLen += int64(len(r.Seq))
    qLen += int64(len(r.Qual))
    for _ ,value := range r.Qual {
        if (value-33) >= 10 && (value-33) < 20 {
	//	fmt.Println(i,value,value-33)
		q10_20 += 1
		
	} 
        if (value-33) >= 20 && (value-33) < 30 {
		q20_30 += 1
	}
	if (value-33) >= 30 {
		q30 += 1
	}
	if (value-33) < 10 {
		q10 += 1
	}
    }
  }
 // fmt.Println("Sample Name:",os.Stdin)
  fmt.Println("Total No of Reads:",n)
  fmt.Println("Total No of Bases:",sLen)
  fmt.Println("Total No of Bases in MB:", (float32(sLen))/1000000)
  fmt.Println("Average Read Length:", (int64(sLen))/(int64(n)))  
 // fmt.Println("Minimum Quality Score:", min)
  fmt.Println("Total No of Bases Q(<10):",q10, (float32(q10))*100/(float32(sLen)))
  fmt.Println("Total No of Bases Q(10-20):",q10_20, (float32(q10_20))*100/(float32(sLen)))
  fmt.Println("Total No of Bases Q(20-30):",q20_30,(float32(q20_30))*100/(float32(sLen)))
  fmt.Println("Total No of Bases Q(>=30):", q30,(float32(q30))*100/(float32(sLen)))
//  fmt.Printf("%v\t%v\t%v\t%v", n, sLen, qLen,q10,q20_30)
}

// Record contains the data from a fasta fastq record
type record struct {
  Name, Seq, Qual string
}

// FqReader holds all the necessary fields that will be use during the processing
// of a fasta fastq file
type FqReader struct {
  Reader          *bufio.Reader
  last, seq, qual []byte // last line processed, temporary seq and qual values
  finished        bool
  rec             record
}

// iterLines iterates over the lines of a reader
func (fq *FqReader) iterLines() ([]byte, bool) {
  line, err := fq.Reader.ReadSlice('\n')
  if err != nil {
    if err == io.EOF {
      return line, true
    } else {
      panic(err)
    }
  }
  return line, false
}

var space = []byte(" ")

func (fq *FqReader) Iter() (record, bool) {
  if fq.finished {
    return fq.rec, fq.finished
  }
  // Read the seq id (fasta or fastq)
  if fq.last == nil {
    for l, done := fq.iterLines(); !done; l, done = fq.iterLines() {
      if l[0] == '>' || l[0] == '@' { // read id
        fq.last = l[0 : len(l)-1]
        break
      }
    }
    if fq.last == nil { // We couldn't find a valid record, no more data in file
      fq.finished = true
      return fq.rec, fq.finished
    }
  }
  fq.rec.Name = string(bytes.SplitN(fq.last, space, 1)[0])
  fq.last = nil

  // Now read the sequence
  fq.seq = fq.seq[:0]
  for l, done := fq.iterLines(); !done; l, done = fq.iterLines() {
    c := l[0]
    if c == '+' || c == '>' || c == '@' {
      fq.last = l[0 : len(l)-1]
      break
    }
    fq.seq = append(fq.seq, l[0:len(l)-1]...)
  }
  fq.rec.Seq = string(fq.seq)

  if fq.last != nil { // There are more lines
    if fq.last[0] != '+' { // fasta record
      return fq.rec, fq.finished
    }
    leng := 0
    fq.qual = fq.qual[:0]
    for l, done := fq.iterLines(); !done; l, done = fq.iterLines() {
      fq.qual = append(fq.qual, l[0:len(l)-1]...)
      leng += len(l)
      if leng >= len(fq.seq) { // we have read enough quality
        fq.last = nil
        fq.rec.Qual = string(fq.qual)
        return fq.rec, fq.finished
      }
    }
    fq.finished = true
    fq.rec.Qual = string(fq.qual)
  }
  return fq.rec, fq.finished // incomplete fastq quality, return what we have
}


