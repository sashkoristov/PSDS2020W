

## Project 1: Object recognition `S3, Lambda Layers, Rekognition, OpenCV` (Java, Python, Node.js)

You'll create a computer vision workflow that intelligently detects dogs, kids (and arbitrary objects) in long streams of video.

#### Motivation

You have a motion-sensing camera installed at your porch. When a motion is detected, the camera starts filming and stores it as separate video file to S3.
At the end of the day you want to know when the dog or kid roamed around too far this day. You can't be bothered to look through the videos manually. 
To complicate things, the camera often activates on birds, neighbors, delivery people.

As most IoT devices, the camera does not do any preprocessing. Therefore the video files are very large. Applying computer vision to the entire videos is infeasible.

Your workflow should efficiently detect at which times of the day your kid or dog appeared. It should be scalable (with number of videos) and cost-efficient (do with as few expensive tasks - Image recognition - as possible).

#### Input 

A folder of videos of that day on S3.

```
└── your-video-bucket
    ├── 1604670378.mp4
    └── 1604670406.mp4
```


#### Rough steps 

* Extract frames (every 1s) from the videos and store them as images to S3. 
* Between each two subsequent images, check if there is a significant delta (something happening). 
* Use AWS Rekognition to detect dogs or kids on those
* Output when either appeared on the porch that day as human sentence (string).



## Project 2: BWA `S3, bwa & samtools, Integrative Genomics Viewer` (Python)


You'll create a real-life workflow to process Escherichia Coli DNA samples, and then investigate whether the patient can be treated with antibiotics.

#### Motivation

You are working in a lab that frequently receives DNA samples of Ecoli bacteria from a hospital.
Ecoli can cause potentially food poisoning, for which the recommended treatment for healthy adults is to just wait it out.
However, in some cases antibiotics may be necessary.

You want to find out if this Ecoli strain could be treated with antibiotics, to then suggest a treatment for the patient.

This has to happen as fast as possible from when the sample arrives. Ecoli DNA is fairly short (4.6 million base pairs). However, samples are usually contaminated (human, bacterial DNA) so processing is hard and has to be parallelized.

Your lab uses a short-read sequencer such as [Illumina MiSeq](https://www.illumina.com/systems/sequencing-platforms/miseq.html) to read ('sequence') the [basepairs](https://en.wikipedia.org/wiki/Base_pair#Examples) of the DNA. 


#### Input


```
└── your-bucket
    ├── NC_000913.3.fasta  
    └── reads
        ├── hipa7_reads_R1.fastq
        └── hipa7_reads_R2.fastq
```

* A [FASTA](https://genome.sph.umich.edu/wiki/FASTA) text file containing the entire DNA of Ecoli ('reference genome')
* A FASTQ text file (FASTA plus likelihood that reads are correct) with paired-end reads ('ABCDE' and 'EDCBA') of your Ecoli sample, obtained from the MiSeq.



#### Rough steps


* Split the reference genome into smaller parts
* Run bwa index 
* Run bwa aln in parallel for each part
* Run bwa sampe in parallel for each part
* Run samtools merge to concat .sam files into one .sam
* Run samtools sort to sort entries
* Run samtools view to convert to .bam (binary sam)
* Run samtools index to make it searchable
* Use IGV to investigate the `hipA` gene for resistance


## Project 3: Prediction of stock prices `S3, Forecast, webhook` (Java, Python, Node.js)

You'll create a real-life workflow to predict the price of stocks you can buy and sell.

#### Motivation
You frequently trade with a wide variety of stocks on various stock exchanges. Every morning, your assistant compiles you a list of 50 stock that have performed interestingly and might be worth a look.

You want to see at a glance how these stocks could perform.
Therefore you program a workflow that takes these names, predicts their prices, and visualises their past and projected price.

This should be done in parallel, per stock. The result should be visualised together.

Furthermore, functions should learn from their RAM usage, so that the workflow becomes more cost-efficient everytime you run it.
This becomes relevant in real-life cases where workflows are run millions of times.


#### Input

```
['AAL', 'ALT', 'NCMI', ... ]
```

* A list of stock ticker symbols traded on some exchange.


#### Rough steps

* Pull commodity prices to S3
* Enter them into AWS Forecast
* Forecast for the coming year for each commodity
* Download the results
* Notify yourself that the prediction is ready


## Prooject 4: Multi-Objective Optimization `Opt4J` (Java)

You will deploy a parallalized optimization algorithm which optimizes a given problem with respect to multiple objectives.

#### Motivation

In this project, you will be working with [Opt4J](http://opt4j.sourceforge.net/), an open-source framework developed for optimization research. Opt4J offers a modular toolbox of diverse optimization components which can be used to implement optimizations with different optimizers, problem representations, evaluation approaches, etc.

In this task, your goal is to (a) study how the optimization of a multi-objective problem via an evolutionary algorithm is implemented within Opt4J, (b) develop a workflow which makes it possible to deploy the evolutionary algorithm as a distributed application, and \(c) evaluate the speedup compared to an execution on your local machine.

#### Input

- Optimizer configuration
- Problem configuration

(A Java project with an optimization problem implemented within Opt4J will be provided as a starting point for this task).

#### Rough steps (Irrelevant for this task)

- Familiarization with the implementation of an evolutionary algorithm within the Opt4J framework. 
- Definition of a parallelization concept: 
	- Dividing the optimization into individual tasks
	- Definition the FC for the evolutionary optimization
- Implementation
	- Implementation of the FC
	- Implementation of the functions


## Project 5: Gate Change Alert `S3, Lambda Rekognition, DB` (Java, Python, Node.js)

You need to create a workflow that performs a series of actions after a gate of a specific flight has changed at an  airport. 

#### Motivation

The workflow reads the information about the flight, the new gate and then loads all available passenger data (from a Database) of that flight. Also, checks the quelength of the security check. Thereafter, for every passenger from that flight that is already at the airport, the workflow reads the gps location and calculates the time to gate.

#### Input 

An image of security check on S3.


#### Rough steps 

- create a database in format you like (in a container, run it in a VM, use some DB service of AWS)
- create a map of the Innsbruck Airport with GPS locations
- distinguish the security and public area
- read the image from S3 and use it as an input for AWS Rekognition


## Parts you can use in addition

### Amazon

```
S3 Glacier
AWS Backup
DynamoDB (NoSQL)
ElastiCache (distributed in-memory DB)
QLDB (immutable transaction log)
Timestream (IoT)
Cloud Map (AWS resource query)
Chatbot (subscribed to SNS)
Polly  (Text-to-speech)
MediaConvert (from * to S3)
Kinesis video streams (from web, sdk, to S3, \*)
Lex (barebones Alexa)
Forecast (time-series)
Rekognition (Object and scene detection, facial, face comparison, text in image)
Comprehend (Natural text processing: entities, key phrases, personal data, language, sentiment)
Translate
Athena (big data query service, basically SQL)
Data exchange (miscellaneous data)
Sumerian (three.js in Browser)
SQS (pub-sub queue)
Greengrass (IoT, local / in-cloud processing)
```

### IBM Lite

```
Watson assistant (chatbot / text to some action)
Visual recognition (scenes, objects...)
Cloudant (JSON database)
Annotator for clinical data
Translate
Object storage
Natural language understanding
Speech to text (streaming transcription)
Text to speech
Tone analyzer (text sentiment)	
```
